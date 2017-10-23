#include "quilting.h"

int main(int argc, char** argv)
{

	Quilting q;

	if (argc == 2)
		q.loadInput(std::string(argv[1]));
	else
		q.loadInput(TEMPO_PATH+"quilting_input.png");

//	std::cout << "Compute random"<< std::endl;
//	q.computeRandom(1024,1024,64);
//	q.saveOutput(TEMPO_PATH+"quilting_random_output.png");

//	std::cout << "Compute fitting"<< std::endl;
//	q.computeFittingTiles(1024,1024,64,8,3);
//	q.saveOutput(TEMPO_PATH+"quilting_fitting_output.png");

	std::cout << "Compute fitting + cutting"<< std::endl;
	q.computeFittingTilesPathCut(1024,1024,64,16,5,4096);
	q.saveOutput(TEMPO_PATH+"quilting_path_output.png");


  return EXIT_SUCCESS;
}

