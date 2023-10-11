pipeline {
    agent any
    stages {
        stage('Clone') {
            steps {
                git branch: 'main', credentialsId: 'cd9bc532-def4-4679-ab1c-87049584e2ea', url: 'https://github.com/quantori/qtr-fingerprint'
            }
        }
        stage('Build') {
            steps {
                sh 'docker build ./python/viz -t qtrshiny'
            }
        }
        stage('Restart') {
            steps {
                sh 'docker stop qtrshiny || true'
                sh 'docker container remove qtrshiny || true'
                sh 'docker run --name qtrshiny --detach -p 8000:8000 --add-host=host.docker.internal:host-gateway -e QTR_API_HOST=http://host.docker.internal qtrshiny'
            }
        }

    }
}