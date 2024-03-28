(part:dev_manual:chap:start_developing_and_contribute:sec:tools:subsec:typos)=
Checking For Typos
=====

When creating a pull request, the Geodynamic World Builder runs a script that checks your code for typos to ensure that new changes are as mistake-free as possible. The tester will output the location of typos after you open a pull request, but to avoid making extra commits you can first check your branch by running the typo test locally. To do this, you will need to install rust and cargo, which can be done following the documentation [here](https://doc.rust-lang.org/cargo/getting-started/installation.html). After installation is complete, run the following command:

`source $HOME/.cargo/env && cargo install typos-cli`

and you can now run `typos` in your `$WORLD_BUILDER_SOURCE_DIR` to check all files within your repository for typos.