#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wall"
#pragma clang diagnostic ignored "-Wextra"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreorder"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#if __GNUG__>5
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#endif
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wreturn-type"
#pragma GCC diagnostic ignored "-Wparentheses"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#elif defined(_MSC_VER)
#pragma warning( push, 0)
#endif


// This code creates pseudo - Wang tiles from input and noise images by inpainting


#include <istream>
#include <fstream>
#include <sstream>
#include "hv/hvPicture.h"
#include <ASTex/special_io.h>
#include <ASTex/utils.h>
#include <string>
#include "hv_astex_io.h"

using namespace hview;

#define LOAD_FILT  // filtered image already exists


void almost_wang_tiles(const std::string& filename_source, const std::string& base_dir, const int nb_clusters)
{
	const int BORDER = 2; // for wang tile generation, apply blending along cuts
	const int WANGBLEND = 1;
	const int NTESTS = 25; // used to generate the tiles
	const int NTESTSI = 20;
	const int maxpatch = 4; // the higher, the more variety in tiles but more time

//	const int CSIZE[nb_clusters] = { 32 }; // used to filter image if no one exists
	int MAX_TP_LEVELS = nb_clusters;
	int MAX_CL = nb_clusters;

	std::string filename_plus_ext = ASTex::IO::remove_path(filename_source);
	std::string name_file = ASTex::IO::remove_ext(filename_plus_ext);

	//	const char *path_in = base_dir.c_str();
	//	const char* path_o = base_dir.c_str();
	//	const char *fname = name_file.c_str();


	int i, j, k;

	hvPictRGB<unsigned char> pioriginal;

	ASTex::load_to_hv(filename_source, pioriginal);

	std::vector<hvPictRGB<unsigned char>> distances(MAX_TP_LEVELS);
	for (k = 0; k < nb_clusters; k++)
	{
		ASTex::load_to_hv(base_dir + name_file + "_mask_bin" + std::to_string(k) + ".png", distances[k]);
	}

	hvPict<unsigned short> origclasses(pioriginal.sizeX(), pioriginal.sizeY(), 0);
	for (i = 0; i < pioriginal.sizeX(); i++) for (j = 0; j < pioriginal.sizeY(); j++)
	{
		int kmax = 0;
		for (k = 1; k < nb_clusters; k++)
		{
			if (distances[kmax].get(i, j).RED() < distances[k].get(i, j).RED())
			{
				kmax = k;
				origclasses.update(i, j, (unsigned short)k);
			}
		}
	}
	// make class colors
	std::vector<hvColRGB<unsigned char> > lcol(MAX_CL);
	std::vector<hvColRGB<double>> coltab(MAX_CL);
	std::vector<int> npixclass(MAX_CL);
	lcol.clear();
	for (i = 0; i < nb_clusters; i++) { coltab[i] = hvColRGB<double>(0.0); npixclass[i] = 0; }
	for (i = 0; i < origclasses.sizeX(); i++) for (j = 0; j < origclasses.sizeY(); j++)
	{
		unsigned short ind = origclasses.get(i, j);
		if (ind >= nb_clusters) hvFatal("incoherent index in makeColorList");
		npixclass[ind]++;
		hvColRGB<unsigned char> col = pioriginal.get(i, j);
		coltab[ind] += (hvColRGB<double>)col;
	}
	for (i = 0; i < nb_clusters; i++)
		lcol.push_back(hvColRGB<unsigned char>((unsigned char)(coltab[i].RED() / (double)npixclass[i]), (unsigned char)(coltab[i].GREEN() / (double)npixclass[i]), (unsigned char)(coltab[i].BLUE() / (double)npixclass[i])));
	// save class pictures
	hvPictRGB<unsigned char> pcl(pioriginal.sizeX(), pioriginal.sizeY(), hvColRGB<unsigned char>(255));
	for (i = 0; i < pioriginal.sizeX(); i++) for (j = 0; j < pioriginal.sizeY(); j++)
	{
		k = origclasses.get(i, j);
		hvColRGB<unsigned char> cc = lcol.at(k);
		pcl.update(i, j, cc);
	}

	// make filtered pictures
	hvPictRGB<unsigned char> fcl, diff, forig;
	hvPict<unsigned char> mask;

	fcl.reset(pioriginal.sizeX(), pioriginal.sizeY(), hvColRGB<unsigned char>(255));
	forig.reset(pioriginal.sizeX(), pioriginal.sizeY(), hvColRGB<unsigned char>(255));
	diff.reset(pioriginal.sizeX(), pioriginal.sizeY(), hvColRGB<unsigned char>(255));
	mask.reset(pioriginal.sizeX(), pioriginal.sizeY(), 0);
	for (i = 0; i < pioriginal.sizeX(); i++) for (j = 0; j < pioriginal.sizeY(); j++)
	{
		k = origclasses.get(i, j);
		int ii, jj;
		for (ii = -BORDER; ii <= BORDER; ii++) for (jj = -BORDER; jj <= BORDER; jj++)
		{
			int px = i + ii, py = j + jj;
			if (px >= 0 && py >= 0 && px < pioriginal.sizeX() && py < pioriginal.sizeY())
			{
				if (origclasses.get(px, py) != k)
				{
					mask.update(i, j, 255);
				}
			}
		}
	}
	if (BORDER > 0)
	{
		mask.gaussianBlur(BORDER, BORDER); mask.gaussianBlur(BORDER, BORDER);
		mask.gaussianBlur(BORDER, BORDER); mask.gaussianBlur(BORDER, BORDER);
	}


	//
	// loading base_dir++name_file+"_Filtered_S.png";
	//
	ASTex::load_to_hv(base_dir + name_file + "_filtered_S.png", forig);

	for (i = 0; i < pioriginal.sizeX(); i++) for (j = 0; j < pioriginal.sizeY(); j++)
	{
		hvColRGB<unsigned char> col = forig.get(i, j);
		col.blend((hvColRGB<double>)pioriginal.get(i, j), col, (double)mask.get(i, j) / 255.0);
		fcl.update(i, j, (hvColRGB<unsigned char>)col);
		hvColRGB<unsigned char> cc; cc.difference(pioriginal.get(i, j), fcl.get(i, j), 127.0, 127.0);
		diff.update(i, j, cc);
	}



	// make the noise images by inpainting
	/////////////////////////////////////////////////////////////
	std::vector<hvPictRGB<unsigned char>> noisei(MAX_CL);
	std::vector<hvPictRGB<unsigned char>> noiseisub(MAX_CL);
	std::vector<hvPictRGB<unsigned char>> noiseimask(MAX_CL);
	int ind;
	for (ind = 0; ind < nb_clusters; ind++)
	{
		noiseimask[ind].reset(pioriginal.sizeX(), pioriginal.sizeY(), hvColRGB<unsigned char>(0));
		noiseisub[ind].reset(pioriginal.sizeX(), pioriginal.sizeY(), hvColRGB<unsigned char>(255));
		hvBitmap map(pioriginal.sizeX(), pioriginal.sizeY(), false);
		for (i = 0; i < origclasses.sizeX(); i++) for (j = 0; j < origclasses.sizeY(); j++)
		{
			if (origclasses.get(i, j) == ind) map.set(i, j, true);
			if (map.get(i, j))
			{
				if (mask.get(i, j) > 0) noiseimask[ind].update(i, j, hvColRGB<unsigned char>(255 - mask.get(i, j)));
				else noiseimask[ind].update(i, j, hvColRGB<unsigned char>(255));
				noiseisub[ind].update(i, j, diff.get(i, j));
			}
		}

		hvColRGB<unsigned char> avgcol = diff.avg(map);
		map.erosion(3, 3); map.erosion(3, 3);
		if (map.count() <= 5)
			noisei[ind].reset(pioriginal.sizeX(), pioriginal.sizeY(), avgcol);
		else
		{
			noisei[ind].reset(pioriginal.sizeX(), pioriginal.sizeY(), hvColRGB<unsigned char>(127));
			noisei[ind].lapped(diff, map, 1.0);
		}
	}


	// make Wang tile of structure
	/////////////////////////////////////////////////////////////

	const double pow = 2.0;
	const double powfilt = 1.0;
	const int blendborder = 3;
	const int blendborderfilt = 3;
	int npatch = 0;

	hvPictRGB<unsigned char> wangpict(pioriginal.sizeX() * 4, pioriginal.sizeY() * 4, hvColRGB<unsigned char>(0));
	hvPictRGB<unsigned char> wangfiltpict(pioriginal.sizeX() * 4, pioriginal.sizeY() * 4, hvColRGB<unsigned char>(0));
	hvPictRGB<unsigned char> wangfiltwpict(pioriginal.sizeX() * 4, pioriginal.sizeY() * 4, hvColRGB<unsigned char>(0));
	std::vector<hvPictRGB<unsigned char>> wangnoiseimask(nb_clusters);
	for (ind = 0; ind < nb_clusters; ind++) wangnoiseimask[ind].reset(pioriginal.sizeX() * 4, pioriginal.sizeY() * 4, hvColRGB<unsigned char>(0));
	hvBitmap wangyet(pioriginal.sizeX() * 4, pioriginal.sizeY() * 4, false);
	for (i = 0; i < pioriginal.sizeX() * 4; i++) for (j = 0; j < pioriginal.sizeY() * 4; j++)
	{
		wangpict.update(i, j, pioriginal.get(i%pioriginal.sizeX(), j%pioriginal.sizeY()));
		wangfiltpict.update(i, j, forig.get(i%pioriginal.sizeX(), j%pioriginal.sizeY()));
		wangfiltwpict.update(i, j, fcl.get(i%pioriginal.sizeX(), j%pioriginal.sizeY()));
		for (ind = 0; ind < nb_clusters; ind++) wangnoiseimask[ind].update(i, j, noiseimask[ind].get(i%pioriginal.sizeX(), j%pioriginal.sizeY()));
	}

	hvVec2<int> fmin, fmax;
	int dimx, dimy;
	dimx = pioriginal.sizeX() / 2;
	dimy = pioriginal.sizeY() / 2;
	hvPict<double> imgHT, imgHB, imgVL, imgVR;
	int *pathHT = new int[dimx];
	int *pathHB = new int[dimx];
	int *pathVL = new int[dimy];
	int *pathVR = new int[dimy];

	hvBitmap fragfeature(dimx, dimy, false);

	// fill from corner
	bool search = true, bordersdone = false, firstcorner = true;
	int cornercc = 0;
	int tilex = 0, tiley = 0, tilecounter = 0;
	int counter = 0;
	do {
		int modulox = 0; int moduloy = 0;
		int locx = pioriginal.sizeX() - 1, locy = pioriginal.sizeY() - 1;
		search = true;

		// first manage to fill the corner and the borders
		if (!bordersdone)
		{
			const int seamlen = 4;

			locx = pioriginal.sizeX() - 1; locy = pioriginal.sizeY() - 1;
			if (!wangyet.get(locx, locy)) search = false;
			if (search)
			{
				locx = pioriginal.sizeX() - 1; locy = pioriginal.sizeY() - 1;
				while (wangyet.get(locx, locy) && locx > pioriginal.sizeX() - 4 && locy > pioriginal.sizeY() - 4) { locx--; locy--; }
				if (!wangyet.get(locx, locy)) { firstcorner = true; search = false; }
			}

			if (search) {
				locx = pioriginal.sizeX() / 2; locy = pioriginal.sizeY() - 1;
				if (!wangyet.get(locx, locy)) search = false;
			}
			if (search) {
				locx = pioriginal.sizeX() - 1; locy = pioriginal.sizeY() / 2;
				if (!wangyet.get(locx, locy)) search = false;
			}
			if (search) {
				locx = pioriginal.sizeX() / 2; locy = 2 * pioriginal.sizeY() - 1;
				if (!wangyet.get(locx, locy)) { search = false; moduloy = 1; }
			}
			if (search) {
				locx = 2 * pioriginal.sizeX() - 1; locy = pioriginal.sizeY() / 2;
				if (!wangyet.get(locx, locy)) { search = false; modulox = 1; }
			}

			if (search)
			{
				locx = pioriginal.sizeX() - 1; locy = pioriginal.sizeY() - 1;
				while (wangyet.get(locx, locy) && locx > 1) locx--;
				int cc = locx - 1; while (cc > 0 && !wangyet.get(cc, locy)) cc--;
				if (locx > 1 && (locx - cc) > seamlen) search = false;
				else if (locx > 1) for (; cc <= locx; cc++) wangyet.set(cc, locy, true);
			}

			if (search)
			{
				locx = pioriginal.sizeX() - 1; locy = pioriginal.sizeY() - 1;
				while (wangyet.get(locx, locy) && locy > 1) locy--;
				int cc = locy - 1; while (!wangyet.get(locx, cc) && cc > 0) cc--;
				if (locy > 1 && (locy - cc) > seamlen) search = false;
				else if (locy > 1) for (; cc <= locy; cc++) wangyet.set(locx, cc, true);
			}

			if (search)
			{
				locx = pioriginal.sizeX() - 1; locy = 2 * pioriginal.sizeY() - 1;
				while (wangyet.get(locx, locy) && locx > 1) locx--;
				int cc = locx - 1; while (!wangyet.get(cc, locy) && cc > 0) cc--;
				if (locx > 1 && (locx - cc) > seamlen) { search = false; moduloy = 1; }
				else if (locx > 1) for (; cc <= locx; cc++) wangyet.set(cc, locy, true);
			}

			if (search)
			{
				locx = 2 * pioriginal.sizeX() - 1; locy = pioriginal.sizeY() - 1;
				while (wangyet.get(locx, locy) && locy > 1) locy--;
				int cc = locy - 1; while (!wangyet.get(locx, cc) && cc > 0) cc--;
				if (locy > 1 && (locy - cc) > seamlen) { search = false; modulox = 1; }
				else if (locy > 1) for (; cc <= locy; cc++) wangyet.set(locx, cc, true);
			}
			while (search && locy > 1)
			{
				locx = 2 * pioriginal.sizeX() - 1; locy = pioriginal.sizeY() - 1;
				while (wangyet.get(locx, locy) && locy > 1) locy--;
				int cc = locy - 1; while (!wangyet.get(locx, cc) && cc > 0) cc--;
				if (locy > 1 && (locy - cc) > seamlen) { search = false; modulox = 1; }
				else if (locy > 1) for (; cc <= locy; cc++) wangyet.set(locx, cc, true);
			}

			if (search) { bordersdone = true; }
		}
		// manage the interior of the tiles

		if (bordersdone)
		{
			if (tilex < 4 && tiley < 4)
			{
				int cc = 0;
				do {
					cc++;
					locx = tilex*pioriginal.sizeX() + pioriginal.sizeX() / 2 + int((2.0*(double)rand() / (double)RAND_MAX - 1.0)*(double)(pioriginal.sizeX()) / 4.0);
					locy = tiley*pioriginal.sizeY() + pioriginal.sizeY() / 2 + int((2.0*(double)rand() / (double)RAND_MAX - 1.0)*(double)(pioriginal.sizeY()) / 4.0);
					modulox = tilex; moduloy = tiley;
				} while (cc < 50 && wangyet.get(locx, locy));
				search = false;
			}
			tilecounter++; if (tilecounter == maxpatch) { tilecounter = 0; tilex++; if (tilex >= 4) { tilex = 0; tiley++; } }
			//search=true;
		}

		if (!search)
		{
			int kk;
			hvVec2<int> minpos, minoffset;
			hvColRGB<double> minerr;
			double minf;
			for (kk = 0; kk < (bordersdone ? NTESTSI : NTESTS); kk++)
			{
				double sum = 0.0;
				fmin = hvVec2<int>(int((double)pioriginal.sizeX()*(double)rand() / (double)RAND_MAX*0.5),
					int((double)pioriginal.sizeY()*(double)rand() / (double)RAND_MAX*0.5));
				// choose position to minimize squared difference
				hvVec2<int> ppos = wangpict.chooseMinSquareDiff(locx - dimx / 2 - dimx / 8, locy - dimy / 2 - dimy / 8, dimx / 4, dimy / 4, fmin.X(), fmin.Y(), dimx, dimy, pioriginal, minerr);

				imgHT.squaredDifference(ppos.X(), ppos.Y() + 3 * dimy / 4, dimx, dimy / 4, wangpict, fmin.X(), fmin.Y() + 3 * dimy / 4, pioriginal);
				imgHT.minPathH(0, dimy / 4, 2, pathHT);
				for (int ind = 0; ind < dimx; ind++) sum += imgHT.get(ind, pathHT[ind]);
				imgHB.squaredDifference(ppos.X(), ppos.Y(), dimx, dimy / 4, wangpict, fmin.X(), fmin.Y(), pioriginal);
				imgHB.minPathH(0, dimy / 4, 2, pathHB);
				for (int ind = 0; ind < dimx; ind++) sum += imgHB.get(ind, pathHB[ind]);
				imgVR.squaredDifference(ppos.X() + 3 * dimx / 4, ppos.Y(), dimx / 4, dimy, wangpict, fmin.X() + 3 * dimx / 4, fmin.Y(), pioriginal);
				imgVR.minPathV(0, dimx / 4, 2, pathVR);
				for (int ind = 0; ind < dimy; ind++) sum += imgVR.get(pathVR[ind], ind);
				imgVL.squaredDifference(ppos.X(), ppos.Y(), dimx / 4, dimy, wangpict, fmin.X(), fmin.Y(), pioriginal);
				imgVL.minPathV(0, dimx / 4, 2, pathVL);
				for (int ind = 0; ind < dimy; ind++) sum += imgVL.get(pathVL[ind], ind);

				if (kk == 0) { minf = sum; minpos = ppos; minoffset = fmin; }
				else if (minf > sum) { minf = sum; minpos = ppos; minoffset = fmin; }
			}
			hvVec2<int> startmin = minoffset;
			// refinement
			for (kk = 0; kk < NTESTSI / 2; kk++)
			{
				double sum = 0.0;
				int fminx = startmin.X() + int((double)pioriginal.sizeX() / 16.0*(1.0 - 2.0*(double)rand() / (double)RAND_MAX));
				int fminy = startmin.Y() + int((double)pioriginal.sizeY() / 16.0*(1.0 - 2.0*(double)rand() / (double)RAND_MAX));
				if (fminx < 0) fminx = 0;
				if (fminx >= pioriginal.sizeX() / 2) fminx = pioriginal.sizeX() / 2 - 1;
				if (fminy < 0) fminy = 0;
				if (fminy >= pioriginal.sizeY() / 2) fminy = pioriginal.sizeY() / 2 - 1;
				fmin = hvVec2<int>(fminx, fminy);
				// choose position to minimize squared difference
				hvVec2<int> ppos = wangpict.chooseMinSquareDiff(locx - dimx / 2 - dimx / 8, locy - dimy / 2 - dimy / 8, dimx / 4, dimy / 4, fmin.X(), fmin.Y(), dimx, dimy, pioriginal, minerr);

				imgHT.squaredDifference(ppos.X(), ppos.Y() + 3 * dimy / 4, dimx, dimy / 4, wangpict, fmin.X(), fmin.Y() + 3 * dimy / 4, pioriginal);
				imgHT.minPathH(0, dimy / 4, 2, pathHT);
				for (int ind = 0; ind < dimx; ind++) sum += imgHT.get(ind, pathHT[ind]);
				imgHB.squaredDifference(ppos.X(), ppos.Y(), dimx, dimy / 4, wangpict, fmin.X(), fmin.Y(), pioriginal);
				imgHB.minPathH(0, dimy / 4, 2, pathHB);
				for (int ind = 0; ind < dimx; ind++) sum += imgHB.get(ind, pathHB[ind]);
				imgVR.squaredDifference(ppos.X() + 3 * dimx / 4, ppos.Y(), dimx / 4, dimy, wangpict, fmin.X() + 3 * dimx / 4, fmin.Y(), pioriginal);
				imgVR.minPathV(0, dimx / 4, 2, pathVR);
				for (int ind = 0; ind < dimy; ind++) sum += imgVR.get(pathVR[ind], ind);
				imgVL.squaredDifference(ppos.X(), ppos.Y(), dimx / 4, dimy, wangpict, fmin.X(), fmin.Y(), pioriginal);
				imgVL.minPathV(0, dimx / 4, 2, pathVL);
				for (int ind = 0; ind < dimy; ind++) sum += imgVL.get(pathVL[ind], ind);

				if (minf > sum) { minf = sum; minpos = ppos; minoffset = fmin; }
			}
			counter++;
			fmin = minoffset;
			imgHT.squaredDifference(minpos.X(), minpos.Y() + 3 * dimy / 4, dimx, dimy / 4, wangpict, fmin.X(), fmin.Y() + 3 * dimy / 4, pioriginal);
			imgHT.minPathH(0, dimy / 4, 2, pathHT);
			imgHB.squaredDifference(minpos.X(), minpos.Y(), dimx, dimy / 4, wangpict, fmin.X(), fmin.Y(), pioriginal);
			imgHB.minPathH(0, dimy / 4, 2, pathHB);
			imgVR.squaredDifference(minpos.X() + 3 * dimx / 4, minpos.Y(), dimx / 4, dimy, wangpict, fmin.X() + 3 * dimx / 4, fmin.Y(), pioriginal);
			imgVR.minPathV(0, dimx / 4, 2, pathVR);
			imgVL.squaredDifference(minpos.X(), minpos.Y(), dimx / 4, dimy, wangpict, fmin.X(), fmin.Y(), pioriginal);
			imgVL.minPathV(0, dimx / 4, 2, pathVL);

			// blend the patch

			bool oncorner = false;
			bool onxa = false, onxb = false, onya = false, onyb = false;
			bool pasted = false;
			fragfeature.clear(false);
			for (i = 0; i < dimx; i++) for (j = 0; j < dimy; j++)
			{
				if (i >= pathVL[j] && i <= pathVR[j] + 3 * dimx / 4 && j >= pathHB[i] && j <= pathHT[i] + 3 * dimy / 4)
				{
					fragfeature.set(i, j, true);
					int vx = minpos.X() + i;  while (vx < 0) vx += pioriginal.sizeX();  vx %= pioriginal.sizeX();
					int vy = minpos.Y() + j;  while (vy < 0) vy += pioriginal.sizeY();  vy %= pioriginal.sizeY();
					if ((vx == 0 || vx == pioriginal.sizeX() - 1) && (vy == 0 || vy == pioriginal.sizeX() - 1)) oncorner = true;
					vx = minpos.X() + i; if (vx < 0) vx += 4 * pioriginal.sizeX();
					if (vx == 0 || vx == pioriginal.sizeX() - 1 || vx == pioriginal.sizeX() || vx == 4 * pioriginal.sizeX() - 1) onxa = true;
					if (vx == 2 * pioriginal.sizeX() - 1 || vx == 2 * pioriginal.sizeX() || vx == 3 * pioriginal.sizeX() - 1 || vx == 3 * pioriginal.sizeX()) onxb = true;
					vy = minpos.Y() + j; if (vy < 0) vy += 4 * pioriginal.sizeY();
					if (vy == 0 || vy == pioriginal.sizeY() - 1 || vy == pioriginal.sizeY() || vy == 4 * pioriginal.sizeY() - 1) onya = true;
					if (vy == 2 * pioriginal.sizeY() - 1 || vy == 2 * pioriginal.sizeY() || vy == 3 * pioriginal.sizeY() - 1 || vy == 3 * pioriginal.sizeY()) onyb = true;
				}
			}
			hvPict<unsigned char> pfeat(fragfeature, blendborder, 255);
			pfeat.gaussianBlur(WANGBLEND, WANGBLEND);
			hvPict<unsigned char> pfeatfilt(fragfeature, blendborderfilt, 255);
			pfeatfilt.gaussianBlur(WANGBLEND, WANGBLEND);

			if (!bordersdone && (oncorner || onxa || onxb || onya || onyb))
			{
				minpos = hvVec2<int>(minpos.X() - modulox*pioriginal.sizeX(), minpos.Y() - moduloy*pioriginal.sizeY());
				for (i = 0; i <= 4; i++) for (j = 0; j <= 4; j++)
				{
					bool doblend = false;
					if (oncorner) { if (firstcorner) doblend = true; }
					else
					{
						if (onxa && !onya && !onyb && (i == 0 || i == 1 || i == 4)) doblend = true;
						if (onxb && !onyb && !onya && (i == 2 || i == 3)) doblend = true;
						if (onya && !onxa && !onxb && (j == 0 || j == 1 || j == 4)) doblend = true;
						if (onyb && !onxb && !onxa && (j == 2 || j == 3)) doblend = true;
					}
					if (doblend)
					{
						wangpict.blendRect(pioriginal.sizeX()*(i - 1) + minpos.X(), pioriginal.sizeY()*(j - 1) + minpos.Y(), fmin.X(), fmin.Y(), dimx, dimy, pioriginal, pfeat, (unsigned char)(255), pow, fragfeature, false);
						wangfiltwpict.blendRect(pioriginal.sizeX()*(i - 1) + minpos.X(), pioriginal.sizeY()*(j - 1) + minpos.Y(), fmin.X(), fmin.Y(), dimx, dimy, fcl, pfeatfilt, (unsigned char)(255), powfilt, fragfeature, false);
						wangfiltpict.blendRect(pioriginal.sizeX()*(i - 1) + minpos.X(), pioriginal.sizeY()*(j - 1) + minpos.Y(), fmin.X(), fmin.Y(), dimx, dimy, forig, pfeatfilt, (unsigned char)(255), powfilt, fragfeature, false);
						for (ind = 0; ind < nb_clusters; ind++)
							wangnoiseimask[ind].blendRect(pioriginal.sizeX()*(i - 1) + minpos.X(), pioriginal.sizeY()*(j - 1) + minpos.Y(), fmin.X(), fmin.Y(), dimx, dimy, noiseimask[ind], pfeatfilt, (unsigned char)(255), powfilt, fragfeature, false);
						wangyet.operatorOr(pioriginal.sizeX()*(i - 1) + minpos.X(), pioriginal.sizeY()*(j - 1) + minpos.Y(), 0, 0, dimx, dimy, fragfeature);
						pasted = true;
						cornercc = 0;
						npatch++;
					}
				}
				if (oncorner && firstcorner) { firstcorner = false; cornercc = 0; }
				if (oncorner && !firstcorner) { if (cornercc >= 5) firstcorner = true; else cornercc++; }
			}
			else if (bordersdone && !(oncorner || onxa || onxb || onya || onyb))
			{

				wangpict.blendRect(minpos.X(), minpos.Y(), fmin.X(), fmin.Y(), dimx, dimy, pioriginal, pfeat, (unsigned char)(255), pow, fragfeature, false);
				wangfiltwpict.blendRect(minpos.X(), minpos.Y(), fmin.X(), fmin.Y(), dimx, dimy, fcl, pfeat, (unsigned char)(255), pow, fragfeature, false);
				wangfiltpict.blendRect(minpos.X(), minpos.Y(), fmin.X(), fmin.Y(), dimx, dimy, forig, pfeat, (unsigned char)(255), pow, fragfeature, false);
				for (ind = 0; ind < nb_clusters; ind++)
					wangnoiseimask[ind].blendRect(minpos.X(), minpos.Y(), fmin.X(), fmin.Y(), dimx, dimy, noiseimask[ind], pfeat, (unsigned char)(255), pow, fragfeature, false);
				wangyet.operatorOr(minpos.X(), minpos.Y(), 0, 0, dimx, dimy, fragfeature);
				pasted = true;
				npatch++;
			}

			if (pasted)
			{
			}
			else if (bordersdone) tilecounter--;
		}
	} while (!search);

	ASTex::save_from_hv(base_dir + name_file + "_wang.png", wangpict);
	ASTex::save_from_hv(base_dir + name_file + "_filtw_wang.png", wangfiltwpict);
	ASTex::save_from_hv(base_dir + name_file + "_filt_wang.png", wangfiltpict);

	for (ind = 0; ind < nb_clusters; ind++)
	{
		ASTex::save_from_hv(base_dir + name_file + "_noise_cl" + std::to_string(ind) + "_wang.png", wangnoiseimask[ind]);
	}

}

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#elif defined(_MSC_VER)
#pragma warning( pop )
#endif
