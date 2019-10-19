# Infotaxis
This is an implementation of infotaxis based on the paper [*‘Infotaxis’ as a strategy for searching without gradients*](https://www.nature.com/articles/nature05464) by Vergassola M, Villermaux E and Shraiman BI. A detailed introduction to infotaxis as a searching strategy is presented in the paper and its supported information.

This implementation includes 2D and 3D versions of infotaxis, allows wind to blow in an aribitrary direction, and enables the searcher to move in eight directions instead of four.

Some examples:

<img src="/finding.png" alt="A searcher is finding the source" width="310"/>
<img src="/found.png" alt="A searcher has found the source" width="310"/>

The moving red dot is the searcher and the big blue dot (whose size varies with source radius) is the source. Brightness in the background indicates the magnitudes of probabilities that the source is at those corresponding locations. White line is the trajectory and dark-red dots on the line are places where the searcher gets patches.
