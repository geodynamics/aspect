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

  ret = system ("cd output-checkpoint_01 ; "
                "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_01.prm "
                " ; "
                " echo 'set Output directory = output1.tmp' "
                " ; "
                " rm -rf output1.tmp ; mkdir output1.tmp "
                ") "
                "| ../../aspect -- >/dev/null ");

  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  ret = system ("cd output-checkpoint_01 ; "
                " rm -rf output2.tmp ; mkdir output2.tmp ; "
                " cp output1.tmp/restart* output2.tmp/");
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;


  std::cout << "* now resuming:" << std::endl;
  ret = system ("cd output-checkpoint_01 ; "
                "(cat " ASPECT_SOURCE_DIR "/tests/checkpoint_01.prm "
                " ; "
                " echo 'set Output directory = output2.tmp' "
                " ; "
                " echo 'set Resume computation = true' "
                ") "
                "| ../../aspect -- >/dev/null");
  if (ret!=0)
    std::cout << "system() returned error " << ret << std::endl;

  std::cout << "* now comparing:" << std::endl;

  ret = system ("cd output-checkpoint_01 ; "
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

  // abort current process:
  exit (0);
  return 42;
}


// run this function by initializing a global variable by it
int i = f();
