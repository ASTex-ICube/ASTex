#include "quilting_all.h"

int main(int argc, char** argv)
{

	QuiltingAll q;

	if (argc < 2)
		return 1;


	q.loadInput(std::string(argv[1]));

	int tw = 64;

	if (argc >= 3)
		tw = atoi(argv[2]);

	if (argc >= 5)
		q.setFirstTile(atoi(argv[3]),atoi(argv[4]));

	if (argc >= 6)
		q.setLevelCompare(atoi(argv[5]));


	q.computeFittingTilesPathCutInterp(tw*4,tw*4,tw,tw/4,tw/8,7,4096);
	q.saveOutput(TEMPO_PATH+std::string(argv[1]));


  return EXIT_SUCCESS;
}

