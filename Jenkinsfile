#!groovy

pipeline {
  agent {
    docker {
        image 'dealii/dealii:v8.5.1-gcc-mpi-fulldepscandi-debugrelease'
	label 'has-docker'
    }
  }

  options {
    timeout(time: 2, unit: 'HOURS') 
  }

  parameters {
    booleanParam(defaultValue: false, description: 'Is the pull request approved for testing?', name: 'TRUST_BUILD')
  }

  stages {
    stage ("info") {
      steps {
        echo "PR: ${env.CHANGE_ID} - ${env.CHANGE_TITLE}"
        echo "CHANGE_AUTHOR_EMAIL: ${env.CHANGE_AUTHOR_EMAIL}"
        echo "CHANGE_AUTHOR: ${env.CHANGE_AUTHOR}"
        echo "CHANGE_AUTHOR_DISPLAY_NAME: ${env.CHANGE_AUTHOR_DISPLAY_NAME}"
        echo "building on node ${env.NODE_NAME}"
      }
    }

    stage ("Check permissions") {
      when {
	allOf {
            environment name: 'TRUST_BUILD', value: 'false' 
            not {branch 'master'}
            not {changeRequest authorEmail: "rene.gassmoeller@mailbox.org"}
            not {changeRequest authorEmail: "heister@clemson.edu"}
            not {changeRequest authorEmail: "bangerth@colostate.edu"}
            not {changeRequest authorEmail: "judannberg@gmail.com"}
            not {changeRequest authorEmail: "ja3170@columbia.edu"}
            not {changeRequest authorEmail: "jbnaliboff@ucdavis.edu"}
            not {changeRequest authorEmail: "menno.fraters@outlook.com"}
            not {changeRequest authorEmail: "a.c.glerum@uu.nl"}
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
        sh 'git diff > changes-astyle.diff'
        archiveArtifacts artifacts: 'changes-astyle.diff', fingerprint: true
        sh '''
	  git diff --exit-code || \
	  { echo "Please check indentation, see artifacts in the top right corner!"; exit 1; }
	  '''
        }
    }

    stage('Build') {
      options {timeout(time: 15, unit: 'MINUTES')}
      steps {
        sh '''
          mkdir -p /home/dealii/build-gcc-fast
          cd /home/dealii/build-gcc-fast
          cmake -G "Ninja" \
	  	-D DEAL_II_CXX_FLAGS='-Werror' \
	  	-D ASPECT_TEST_GENERATOR=Ninja \
		-D ASPECT_USE_PETSC=OFF \
		-D ASPECT_RUN_ALL_TESTS=ON \
		-D ASPECT_PRECOMPILE_HEADERS=ON \
		$WORKSPACE/
          ninja
        '''
      } 
    }

    stage('Run tests') {
      options {timeout(time: 90, unit: 'MINUTES')}
      steps {
        sh '''
          rm -f /home/dealii/build-gcc-fast/FAILED
          cd /home/dealii/build-gcc-fast/tests
          echo "Prebuilding tests..."
          ninja -k 0 tests >/dev/null || { touch /home/dealii/build-gcc-fast/FAILED; }
          cd ..
          echo "Checking for test failures..."
          ctest --output-on-failure -j4 || { echo "At least one test FAILED"; }

          echo "Generating reference output..."
          ninja generate_reference_output
        '''
        sh 'git diff tests > changes-test-results.diff'
        archiveArtifacts artifacts: 'changes-test-results.diff', fingerprint: true
        sh 'if [ -f /home/dealii/build-gcc-fast/FAILED ]; then exit 1; fi'
        sh 'git diff --exit-code --name-only'
      }
    }
  }
}
