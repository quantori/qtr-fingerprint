pipeline {
  agent {
    docker {
      image 'rikorose/gcc-cmake:latest'
      args '-v /home/centos/SFO:/home/centos/workspace/SFO/slow-tests/SFO'
      label 'qtr'
    }
  }
  stages {
    stage('build') {
      steps {
        sh '''
          pwd
          ls -lha
          ls -lha ./SFO/
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
          ./bin/tests --gtest_output="xml:./report.xml" --big_data_dir_path=/home/centos/workspace/SFO/slow-tests/SFO --gtest_filter='SlowTest*'
        '''
      }
    }
  }
  post {
      always {
        junit 'cpp/build/report.xml'
      }
  } 
}
