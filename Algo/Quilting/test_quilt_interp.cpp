#include "quilting_interp.h"

int main(int argc, char** argv)
{

	QuiltingInterp q;

	if (argc == 2)
		q.loadInput(std::string(argv[1]));
	else
		q.loadInput(TEMPO_PATH+"quilting_input.png");

	q.computeFittingTilesPathCutInterp(1024,1024,128,16,7,5,4096);
	q.saveOutput(TEMPO_PATH+"quilting_path_output.png");


  return EXIT_SUCCESS;
}

