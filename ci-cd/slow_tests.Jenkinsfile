pipeline {
  agent {
    docker {
      image 'rikorose/gcc-cmake:latest'
      args '-v /home/centos/SFO:/home/centos/workspace/SFO/slow-tests/SFO'
      label 'qtr'
    }
  }
  options {
    buildDiscarder(logRotator(numToKeepStr: '60', daysToKeepStr: '60',  artifactNumToKeepStr: '60', artifactDaysToKeepStr: '60'))
  }
  stages {
    stage('build') {
      steps {
        sh '''
          cd cpp
          mkdir -p build && cd build
          cmake ..
          cmake --build . -j10
        '''
      }
    }
    stage('test') {
      steps {
        sh '''
          cd cpp/build
          export GLOG_logtostdout=true
          ./bin/tests --gtest_output="xml:./report.xml" --big_data_dir_path=/home/centos/workspace/SFO/slow-tests/SFO --gtest_filter='SlowTest*'
        '''
      }
    }
  }
  post {
      always {
        junit 'cpp/build/report.xml'
        cleanWs()
      }
  } 
}
