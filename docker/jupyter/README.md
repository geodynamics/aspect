Jupyter Notebook docker container
---------------------------------

to make the image run

    make

To run a local jupyter notebook server run

    make test

and browse to http://localhost:8888/

Alternatively, you can run

    docker run -d -p 8888:8888 --name tmpnb-aspect-jupyter tjhei/aspect-jupyter start-notebook.sh --NotebookApp.token=''
