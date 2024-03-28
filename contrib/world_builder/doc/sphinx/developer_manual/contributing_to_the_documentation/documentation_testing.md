(part:dev_manual:chap:contribute_to_doc:sec:doc_testing)=
Documentation testing
=====================

Making sure that the examples in the documentation are working is important to make sure users can rely on the documentation. For that reason, all files in `doc/sphinx/_static/gwb_input_files/` which both have a `.wb` and `.grid` file, are automatically tested in every pull request. They are only tested to see if any error occurs, not if the output is correct or has changed.