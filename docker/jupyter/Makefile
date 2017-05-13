TAG=tjhei/aspect-jupyter

img:
	docker build -t ${TAG} .

push: img
	docker push ${TAG}

test: kill
	docker run -d -p 8888:8888 --name tmpnb-aspect-jupyter ${TAG} start-notebook.sh --NotebookApp.token=''
	@echo visit http://localhost:8888 now and type "make kill" to remove

kill:
	docker rm -f tmpnb-aspect-jupyter || true

.PHONY: img test kill
