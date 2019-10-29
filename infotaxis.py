"""
Author: Ziyi Gong
"""

import numpy as np
from scipy.special import k0
from itertools import product


class Infotaxis:
    def __init__(self, d, r, a, v, winddir, tau, dt, dim, pos, src_pos, src_radius, window):
        """
		Parameters
		----------
		d: float
		    Diffusivity (m**2/s)
		r: float
		    Emission rate (particles/s)
		a: float
		    Size of the searcher (m)
		v: float
		    Wind speed (m/s)
		winddir: float
		    Wind direction, between 0 (East) and 2*PI, counterclockwise
		tau: float
		    Particle lifetime (s)
		dim: int
		    Dimension, either 2 or 3
		pos: array-like
		    The position of the searcher
		src_pos: array-like
		    The position of the source, i.e. r0 in the paper
		src_radius: int
		    The radius of the source
		window: int
		    Window size in which the searcher runs
		"""

        self._d = d
        self._r = r
        self._a = a
        self._tau = tau
        self._dt = dt
        self._dim = dim
        self._window = window

        # Lambda in the concentration and hit rate equations
        self._lam = np.sqrt(d * tau / (1 + v ** 2 * tau / (4 * d)))

        # The scalar in the concentration equation
        # r / (2 * PI * d) for dim == 2 and r / (4 * PI * d) for dim == 3
        self._sc = r / (2 * np.pi * d * (dim - 1))

        # The scalar in the hit rate equation
        # 2 * PI * d / ln(lambda / a) for dim == 2 and 4 * PI * d for dim == 3
        if dim == 2:
            self._sr = 2 * np.pi * d / np.log(self._lam / a)
        elif dim == 3:
            self._sr = 4 * np.pi * d
        else:
            raise ValueError('Unsupported dimension')

        self._v = np.array([-np.sin(winddir), np.cos(winddir)]) * v  # Vectorize
        self._v[np.abs(self._v) < 1e-5] = 0

        self._pos = self.xy2ij(pos)

        self._src_pos = self.xy2ij(src_pos)
        self._src_radius = src_radius

        # Distance-related metrics
        self._ticks = np.arange(window, dtype=np.float64)
        self._diffi = None
        self._diffj = None
        self._dists = None

        # The probability of source locating in each position
        self._prob_map = np.full((window, window), 1 / window ** 2)

        self._s = self.entropy()

        # Auxiliary arrays
        self._delta_s_e = np.zeros((3, 3))
        self._hit_rates = np.zeros((3, 3, self._window, self._window))

    def xy2ij(self, pos):
        ij = np.zeros(2, dtype=int)
        ij[0] = self._window - 1 - pos[1]
        ij[1] = pos[0]
        return ij

    def ij2xy(self, indices):
        xy = np.zeros(2, dtype=int)
        xy[0] = indices[1]
        xy[1] = self._window - 1 - indices[0]
        return xy

    def poisson(self, hit_rates):
        """
        Generate a poisson random number and reformat it to 0 (no hit) or 1 (some hits).

		Parameters
		---------
		hit_rates: numpy.ndarray
		    The hit rates for all possible source locations

		Returns
		-------
		hits: int
            0 (no hit) or 1 (some hits)
		"""
        return int(np.random.poisson(hit_rates * self._dt) > 0)

    def entropy(self, probs=None):
        if probs is None:
            probs = self._prob_map
        probs[probs == 0] = 1

        return -np.sum(probs * np.log(probs))

    def adjacent(self, a):
        if a == 0:
            return 1, 2
        elif a == self._window - 1:
            return 0, 1
        else:
            assert (0 < a < self._window)
            return 0, 1, 2

    @property
    def cur(self):
        return self.ij2xy(self._pos), self._a

    @property
    def src(self):
        return self.ij2xy(self._src_pos), self._src_radius

    @property
    def prob_map(self):
        return self._prob_map

    @property
    def wind(self):
        return self._v[1], -self._v[0]

    def concentration_map(self):
        """
        Returns
        -------
        concentrations: numpy.ndarray
		    The concentration corresponding to all possible source locations when the searcher
		    is at this position
		"""
        exp_part = np.exp(np.add.outer(self._diffi * self._v[0], self._diffj * self._v[1])
                          / (2 * self._d))

        self._dists[self._dists == 0] += 0.1

        if self._dim == 2:
            rest = k0(self._dists / self._lam)
        else:
            rest = np.exp(-self._dists / self._lam) / self._dists

        rest[self._dists < self._src_radius] = np.inf

        return self._sc * exp_part * rest

    def rate_map(self, concentration):
        """
        Returns
        -------
        rates: numpy.ndarray
            The concentration corresponding to all possible source locations when the searcher
            is at this position
        """
        hit_rates = self._sr * concentration
        return hit_rates

    def update_posterior(self, hit_rates, num_hits, inplace=False):
        """
		Update the posterior probability distribution 
		P(src_pos|trajectory(t))
			= normalize(P(src_pos | trajectory(t-1)) * P(hits(self._pos(t)) | src_pos))
		Parameters
		----------
		hit_rates: numpy.ndarray
		    The hit rates for all possible source locations
		num_hits: int
		    The number of hits encountered at t. It should between 0 (no hit) and 1 (hits)
		"""
        p_h = np.exp(-hit_rates * self._dt)
        if num_hits == 1:
            p_h = 1 - p_h

        if inplace:
            self._prob_map *= p_h
            self._prob_map /= self._prob_map.sum()

        else:
            new_probs = self._prob_map * p_h
            return new_probs / new_probs.sum()

    def P_miss(self, hit_rates):
        """
        P_miss = sum(P(h(self._pos)==0|src_pos) * P(src_pos|trajectory(t)))

        Parameters
        ----------
        hit_rates: numpy.ndarray
		    The hit rates for all possible source locations

        Returns
        -------
        p_miss: float
            The probability of getting no hit
        """

        return np.sum(np.exp(- self._dt * hit_rates) * self._prob_map)

    def step(self):
        """
        Take a step. The step can be one of the 9 pixels in a 3-by-3 grid centered at the
        searcher's current discretized location.

        Returns
        -------
        pos: numpy.ndarray
            The next position
        hit: bool
            Whether the searcher gets patches or not
        """
        # Check if source has been found
        if np.linalg.norm(self._pos - self._src_pos) < self._src_radius:
            return -1

        # Clear
        self._delta_s_e[:, :] = 0
        self._hit_rates[:, :, :, :] = np.inf

        # Get the possible neighbor indices
        xs, ys = self.adjacent(self._pos[0]), self.adjacent(self._pos[1])

        for i, j in product(xs, ys):  # For each neighbor
            # Calculate distances to all locations
            self._diffi = self._pos[0] + (i - 1) - self._ticks
            self._diffj = self._pos[1] + (j - 1) - self._ticks
            self._dists = np.add.outer(self._diffi ** 2, self._diffj ** 2)
            self._dists = np.sqrt(self._dists)

            # Get the concentration map
            concentrations = self.concentration_map()
            self._hit_rates[i, j] = self.rate_map(concentrations)

            # Get probability of finding source
            p_found = self._prob_map[self._dists < self._src_radius].sum()
            p_not_found = 1 - p_found

            # If it does not find the source:

            # Get the probability of missing/getting a hit
            p_miss = self.P_miss(self._hit_rates[i, j])
            p_hit = 1 - p_miss

            s_hit = self.entropy(self.update_posterior(self._hit_rates[i, j], 1))
            s_miss = self.entropy(self.update_posterior(self._hit_rates[i, j], 0))

            # get expected entropy decrease given source not found
            delta_s_not_found = p_hit * (s_hit - self._s) + p_miss * (s_miss - self._s)

            # compute total expected entropy decrease
            self._delta_s_e[i, j] = p_found * -self._s + p_not_found * delta_s_not_found

        # Choose move that decreases entropy the most
        argmin_1d = np.argmin(self._delta_s_e)
        argmin = np.zeros(2, dtype=int)
        argmin[0] = argmin_1d // 3
        argmin[1] = argmin_1d % 3
        self._pos[0] += argmin[0] - 1
        self._pos[1] += argmin[1] - 1

        # Check if source has been found
        if np.linalg.norm(self._pos - self._src_pos) < self._src_radius:
            self._prob_map[:, :] = 0
            return -1

        # Get the hits at the current position given hit rates
        hit = self.poisson(self._hit_rates[argmin[0], argmin[1],
                                           self._src_pos[0], self._src_pos[1]])

        # Update source posterior
        self.update_posterior(self._hit_rates[argmin[0], argmin[1]], hit, inplace=True)

        return self.ij2xy(self._pos), bool(hit)
