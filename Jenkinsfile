#!groovy

// load library https://github.com/tjhei/jenkins-stuff to provide
// killold.killOldBuilds() function:
@Library('tjhei') _

pipeline {
  agent none

  options {
    timeout(time: 3, unit: 'HOURS')
  }

  stages
  {
    stage("abort old")
    {
      agent none
      steps
      {
        githubNotify context: 'Jenkins', description: 'initializing...',  status: 'PENDING'
        // kill older builds in this PR:
        script { killold.killOldBuilds() }
      }
    }

    stage("main")
    {
      agent
      {
        // first build a docker image for running the tests from the repo
        docker {
          image 'geodynamics/aspect-tester:focal-dealii-9.5-v3'
          // We mount /repos into the docker image. This allows us to cache
          // the git repo by setting "advanced clone behaviors". If the
          // directory does not exist, this will be ignored.
          args '-v /repos:/repos:ro'
        }
      }

      post { cleanup { cleanWs() } }

      stages {
        stage ("Print Info") {
          steps {
            echo "PR: ${env.CHANGE_ID} - ${env.CHANGE_TITLE}"
            echo "CHANGE_AUTHOR_EMAIL: ${env.CHANGE_AUTHOR_EMAIL}"
            echo "CHANGE_AUTHOR: ${env.CHANGE_AUTHOR}"
            echo "CHANGE_AUTHOR_DISPLAY_NAME: ${env.CHANGE_AUTHOR_DISPLAY_NAME}"
            echo "building on node ${env.NODE_NAME}"
          }
        }

        stage ("Check Permissions") {
          when {
            // check for "ready to test" if it is a PR and not by one of the people listed
            allOf {
              changeRequest()
              not {changeRequest authorEmail: "rene.gassmoeller@mailbox.org"}
              not {changeRequest authorEmail: "timo.heister@gmail.com"}
              not {changeRequest authorEmail: "bangerth@colostate.edu"}
              not {changeRequest authorEmail: "judannberg@gmail.com"}
              not {changeRequest authorEmail: "ja3170@columbia.edu"}
              not {changeRequest authorEmail: "jbnaliboff@ucdavis.edu"}
              not {changeRequest authorEmail: "menno.fraters@tutanota.com"}
              not {changeRequest authorEmail: "a.c.glerum@uu.nl"}
              not {changeRequest authorEmail: "myhill.bob@gmail.com"}
              not {changeRequest authorEmail: "ljhwang@ucdavis.edu"}
            }
          }
          steps {
            // For /rebuild to work you need to:
            // 1) select "issue comment" to be delivered in the github webhook setting
            // 2) install "GitHub PR Comment Build Plugin" on Jenkins
            // 3) in project settings select "add property" "Trigger build on pr comment" with
            //    the phrase ".*/rebuild.*" (without quotes)
            sh '''
            wget -q -O - https://api.github.com/repos/geodynamics/aspect/issues/${CHANGE_ID}/labels | grep 'ready to test' || \
            { echo "This commit will only be tested when it has the label 'ready to test'. Trigger a rebuild by adding a comment that contains '/rebuild'..."; exit 1; }
            '''
          }
          post {
            failure
            {
              githubNotify context: 'Jenkins', description: 'need ready to test label and /rebuild',  status: 'PENDING'
              script
              {
                currentBuild.result='NOT_BUILT'
              }
            }
          }
        }

        stage('Check Indentation') {
          steps {
            sh './contrib/utilities/indent'
            sh 'git diff > changes-astyle.diff'
            archiveArtifacts artifacts: 'changes-astyle.diff', fingerprint: true
            sh '''
            git diff --exit-code || \
            { echo "Please check indentation, see artifacts in the top right corner!"; exit 1; }
            '''
          }
          post
          {
            failure
            {
              githubNotify context: 'Jenkins', description: 'indentation failed',  status: 'PENDING'
            }
          }
        }

        stage('Build') {
          options {
            timeout(time: 30, unit: 'MINUTES')
          }
          steps {
            sh '''
            # running cmake...
            mkdir build
            cd build

            cmake \
            -G 'Ninja' \
            -D CMAKE_CXX_FLAGS='-Werror' \
            -D CMAKE_BUILD_TYPE='DebugRelease' \
            -D ASPECT_ADDITIONAL_CXX_FLAGS='-O3' \
            -D ASPECT_TEST_GENERATOR='Ninja' \
            -D ASPECT_PRECOMPILE_HEADERS=ON \
            -D ASPECT_UNITY_BUILD=ON \
            -D ASPECT_WITH_NETCDF=ON \
            -D ASPECT_RUN_ALL_TESTS='ON' \
            -D ASPECT_INSTALL_EXAMPLES='ON' \
            ..
            '''

            sh '''
            # compiling...
            cd build
            ninja
            '''
          }
          post {
            always
            {
              archiveArtifacts artifacts: 'build/detailed.log', fingerprint: true
            }
            failure
            {
              githubNotify context: 'Jenkins', description: 'build failed',  status: 'FAILURE'
            }
          }
        }

        stage('Build Documentation') {
          steps {
            sh 'cd doc && ./update_parameters.sh ./build/aspect'
          }
          post {
            failure
            {
              githubNotify context: 'Jenkins', description: 'documentation failed',  status: 'FAILURE'
            }
          }
        }

        stage('cookbooks') {
          options {
            timeout(time: 20, unit: 'MINUTES')
          }
          steps {
            sh '''
            export BUILDDIR=`pwd`/build
            export NP=`grep -c ^processor /proc/cpuinfo`
            cd cookbooks && make -f check.mk CHECK=--validate BUILD=$BUILDDIR -j $NP
            cd ..
            cd benchmarks && make -f check.mk CHECK=--validate BUILD=$BUILDDIR -j $NP
            '''
          }
          post {
            failure
            {
              githubNotify context: 'Jenkins', description: 'cookbooks failed',  status: 'FAILURE'
            }
          }
        }

        stage('Test') {
          options {
            timeout(time: 150, unit: 'MINUTES')
          }
          steps {
            // Let ninja prebuild the test libraries and run
            // the tests to create the output files in parallel. We
            // want this to always succeed, because it does not generate
            // useful output (we do this further down using 'ctest', however
            // ctest can not run ninja in parallel, so this is the
            // most efficient way to build the tests).
            sh '''
            # prebuilding tests...
            cd build/tests
            ninja -k 0 tests || true
            '''

            // Output the test results using ctest. Since
            // the tests were prebuild in the previous shell
            // command, this will be fast although it is not
            // running in parallel. We use the ninja test generator, which
            // does not support running ctest with -j.
            sh '''
            # generating test results...
            cd build
            ctest \
            --no-compress-output \
            --test-action Test \
            --output-on-failure
            '''
          }
          post {
            always {
              // Generate the 'Tests' output page in Jenkins
              xunit testTimeMargin: '3000',
              thresholdMode: 1,
              thresholds: [failed(), skipped()],
              tools: [CTest(pattern: 'build/Testing/**/*.xml')]

              // Update the reference test output with the new test results
              sh '''
              # generating reference output...
              cd build
              ninja generate_reference_output
              '''

              // Generate the 'Artifacts' diff-file that can be
              // used to update the test results
              sh 'git diff tests > changes-test-results.diff'
              archiveArtifacts artifacts: 'changes-test-results.diff'

              githubNotify context: 'Jenkins', description: 'OK',  status: 'SUCCESS'

            }
            failure
            {
              githubNotify context: 'Jenkins', description: 'tests failed',  status: 'FAILURE'
            }
          }
        }
      }
    }
  }
}
