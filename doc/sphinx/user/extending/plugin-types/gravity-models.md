# Gravity models

The gravity model is responsible for describing the magnitude and direction of
the gravity vector at each point inside the domain. To implement a new gravity
model, you need to overload the [aspect::GravityModel::Interface][20] class
and use the `ASPECT_REGISTER_GRAVITY_MODEL` macro to register your new class.
The implementation of the new class should be in namespace
`aspect::GravityModel`.

Specifically, your new class needs to implement the following basic interface:

```{code-block} c++
template <int dim>
    class aspect::GravityModel::Interface
    {
      public:
        virtual
        Tensor<1,dim>
        gravity_vector (const Point<dim> &position) const = 0;

        virtual
        void
        update ();

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);
    };
```

The kind of information these functions need to provide is discussed in the
documentation of this interface class at
[aspect::GravityModel::Interface][20]. The first needs to return a gravity
vector at a given position, whereas the second is called at the beginning of
each time step, for example to allow a model to update itself based on the
current time or the solution of the previous time step. The purpose of the
last two functions has been discussed in the general overview of plugins
above.
