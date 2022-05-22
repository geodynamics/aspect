# Parameter file

You can create the parameter file in the same way you would do it for any
other model. Beyond that, make sure to start the file with a comment that
explains what this cookbook is about in a few sentences. After that, you will
list all of the input parameters. In general, it makes sense to begin with the
ones that are most important for the setup you want to show, and otherwise to
group parameters and sections that are related to each other (like all
boundary conditions or all initial conditions). To make the input file easy to
understand for other users, it is a good practice to add a short comment to
each section or important parameter used in the file, explaining what this
input option accomplishes and why it is needed for the model setup.

Once you have finalized your input file, you can put it into the [cookbooks](https://github.com/geodynamics/aspect/tree/main/cookbooks)
folder.
