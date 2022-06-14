pipeline {
  agent {
    docker {
      image 'rikorose/gcc-cmake:latest'
      label 'qtr'
    }
  }
  stages {
    stage('build') {
      steps {
        sh '''
          cd cpp
          mkdir build && cd build
          cmake ..
          cmake --build . -j10
        '''
      }
    }
    stage('test') {
      steps {
        sh '''
          cd cpp/build
          ./bin/tests --gtest_output="xml:./report.xml"
        '''
        publishHTML target: [
            allowMissing: false,
            alwaysLinkToLastBuild: false,
            keepAll: true,
            reportDir: 'bin',
            reportFiles: 'report.html',
            reportName: 'Slow test report'
          ]
      }
    }
  }
  post {
    cleanup {
      deleteDir()
    }
  }
}
