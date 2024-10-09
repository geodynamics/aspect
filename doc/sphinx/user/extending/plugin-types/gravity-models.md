# Gravity models

The gravity model is responsible for describing the magnitude and direction of
the gravity vector at each point inside the domain. To implement a new gravity
model, you need to overload the [aspect::GravityModel::Interface](https://aspect.geodynamics.org/doc/doxygen/namespaceaspect_1_1GravityModel.html)
class and use the `ASPECT_REGISTER_GRAVITY_MODEL` macro to register your new class.
The implementation of the new class should be in namespace
`aspect::GravityModel`.

The functions you need to overload are extensively
discussed in the documentation of this interface class at
[aspect::GravityModel::Interface](https://aspect.geodynamics.org/doc/doxygen/namespaceaspect_1_1GravityModel.html).
