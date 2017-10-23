#include "quilting_interp_cont.h"

int main(int argc, char** argv)
{

	QuiltingInterpCont q;

	if (argc == 2)
		q.loadInput(std::string(argv[1]));
	else
		q.loadInput(TEMPO_PATH+"quilting_input.png");

//	q.setFirstTile(600,90);

	q.computeFittingTilesPathCutInterp(512,512,64,16,7,5,4096);
	q.saveOutput(TEMPO_PATH+"quilting_path_output.png");


  return EXIT_SUCCESS;
}

