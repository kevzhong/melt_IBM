
# Home

## Lagrangian melting front-tracking IBM documentation

This documentation provides a few technical details on the implementation of the `melt_IBM` code.

Most of the technical details on the method are prepared in the following article:

* Zhong, K., Howland, C. J., Lohse, D., and Verzicco, R. (2025). A front-tracking immersed-boundary framework for simulating Lagrangian melting problems. J. Comput. Phys., page 113762.

Some of the details in the article will be repeated here in this documentation for completeness. But the main purpose of this documentation is to highlight some more specific details not mentioned in the article (mostly specific to implementation) and to walkthrough some of the code.

A to-do list of documentation I intend to add:

- Quick-start guide, explanation of input file parameters
- Configuring the HIT forcing, particularly the stochastic scheme of Eswaran & Pope (1988)
- Initial condition configuration, particularly for restarting simulations from a pre-cursor HIT simulation
- Explanation of the .gts file convention for triangulated geometries
- Cell-tagging details: ray-tagging, GenSDF (Patil et al. 2025)
- Step-by-step examples on the remeshing algorithm at work
- Details on the rigid body routines: quaternion usage, hydrodynamic load calculations via MLS interpolation