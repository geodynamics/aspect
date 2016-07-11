#include <aspect/simulator.h>

/*
 * Launch the following function when this plugin is created. Launch ASPECT
 * twice to test checkpoint/resume and then abort the outer ASPECT run.
 */
int f()
{
  std::cout << "* starting from beginning:" << std::endl;

  // call ASPECT with "--" and pipe an existing input file into it.
  int ret;

  ret = system ("cd output-checkpoint_02_petsc ; "
                "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_02_petsc.prm "
                " ; "
                " echo 'set Output directory = output1.tmp' "
                " ; "
                " rm -rf output1.tmp ; mkdir output1.tmp "
                ") "
                "| ../../aspect -- >/dev/null ");

  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  ret = system ("cd output-checkpoint_02_petsc ; "
                " rm -rf output2.tmp ; cp -r output1.tmp output2.tmp ;"
                " rm -f output1.tmp/log.txt; ");
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;


  std::cout << "* now resuming:" << std::endl;
  ret = system ("cd output-checkpoint_02_petsc ; "
                "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_02.prm "
                " ; "
                " echo 'set Output directory = output1.tmp' "
                " ; "
                " echo 'set Resume computation = true' "
                " ; "
                " echo 'subsection Checkpointing' "
                " ; "
                " echo 'set Restart from rotating checkpoint slot = 0' "
                " ; "
                " echo 'end' "
                ") "
                "| ../../aspect -- >/dev/null");

  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  std::cout << "* now comparing:" << std::endl;

  ret = system ("cd output-checkpoint_02_petsc ; "
                "cp output1.tmp/depth_average.gnuplot depth_average.gnuplot1;"
                "cp output2.tmp/depth_average.gnuplot depth_average.gnuplot2;"
                "cp output1.tmp/solution-00009.0000.gnuplot solution-00009.0000.gnuplot1;"
                "cp output2.tmp/solution-00009.0000.gnuplot solution-00009.0000.gnuplot2;"
                "cp output1.tmp/log.txt log.txt1;"
                "cp output2.tmp/log.txt log.txt2;"
                "cp output1.tmp/statistics statistics1;"
                "cp output2.tmp/statistics statistics2;"
                "");
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  ret = system ("cd output-checkpoint_02_petsc ; "
                "diff -c output?.tmp/rotating-checkpoint.resume-1.z;"
                "diff -c output?.tmp/rotating-checkpoint.mesh-1;"
                "");
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  // abort current process:
  exit (0);
  return 42;
}


// run this function by initializing a global variable by it
int i = f();
