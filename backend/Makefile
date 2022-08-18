
init:
	cp .sample.env .env
	cp docker-compose.local.sample.yml docker-compose.local.yml

up:
	docker-compose --file docker-compose.local.yml up --build

down:
	docker-compose --file docker-compose.local.yml down

lint:
	isort ./src --settings-path=.isort.cfg

test:
	python -m pytest
