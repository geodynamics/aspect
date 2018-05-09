#!groovy

pipeline {
  agent {
    docker {
        image 'dealii/dealii:v8.5.1-gcc-mpi-fulldepscandi-debugrelease'
	label 'has-docker'
    }
  }

  parameters {
    booleanParam(defaultValue: false, description: 'Is the pull request approved for testing?', name: 'TRUST_BUILD')
  }

  stages {
    stage ("Check execution") {
      when {
	allOf {
            environment name: 'TRUST_BUILD', value: 'false' 
            not {branch 'master'}
            not {changeRequest authorEmail: "rene.gassmoeller@mailbox.org"}
            not {changeRequest authorEmail: "timo.heister@gmail.com"}
            not {changeRequest authorEmail: "bangerth@colostate.edu"}
            not {changeRequest authorEmail: "judannberg@gmail.com"}
            not {changeRequest authorEmail: "ja3170@columbia.edu"}
            not {changeRequest authorEmail: "john.naliboff@gmail.com"}
            not {changeRequest authorEmail: "menno.fraters@outlook.com"}
            not {changeRequest authorEmail: "acglerum@gfz-potsdam.de"}
	    }
      }
      steps {
	  echo "Please ask an admin to rerun jenkins with TRUST_BUILD=true"
	    sh "exit 1"
      }
    }

    stage('Check indentation') {
      steps {
        sh './doc/indent'
        sh 'git diff | tee changes-astyle.diff'
        archiveArtifacts artifacts: 'changes-astyle.diff', fingerprint: true
        sh 'git diff --exit-code --name-only'
        }
    }

    stage('Build') {
      steps {
        sh '''
          env
          mkdir -p /home/dealii/build-gcc-fast
          cd /home/dealii/build-gcc-fast
          cmake -G "Ninja" gcc -D ASPECT_TEST_GENERATOR=Ninja -D ASPECT_USE_PETSC=OFF -D ASPECT_RUN_ALL_TESTS=ON -D ASPECT_PRECOMPILE_HEADERS=ON $WORKSPACE/
          ninja
        '''
      } 
    }

    stage('Run tests') {
      steps {
        sh '''
          cd /home/dealii/build-gcc-fast/tests
          echo "prebuilding tests..."
          ninja -k 0 tests >/dev/null 2>&1
          cd ..
          echo "+ ctest --output-on-failure -j $NPROC"
          ctest --output-on-failure -j $NPROC || { echo "test FAILED"; }

          echo "+ ninja generate_reference_output"
          ninja generate_reference_output
          echo "ok"
        '''
        archiveArtifacts artifacts: '/home/dealii/build-gcc-fast/changes.diff', fingerprint: true
        sh 'git diff --exit-code --name-only'
      }
    }
  }
}
