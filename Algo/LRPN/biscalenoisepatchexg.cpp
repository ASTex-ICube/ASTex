/*******************************************************************************
* ASTex:                                                                       *
* Copyright (C) IGG Group, ICube, University of Strasbourg, France             *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: https://astex-icube.github.io                                      *
* Contact information: astex@icube.unistra.fr                                  *
*                                                                              *
*******************************************************************************/



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

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#include "hv/Fragmenter.h"
#include "hv/hvProcWPhase.h"
#include <ASTex/special_io.h>
#include <ASTex/utils.h>
#include <string>

#include "hv/hvIntersection3.h"
#include "hv/hvPerspective.h"
#include "hv/hvPlan3.h"
#include "hv_astex_io.h"

using namespace hview;
using namespace Fragmentation;



void Poisson(int sx, int sy, int NSEEDS, int seedx[], int seedy[])
{
	for (int i = 0; i < NSEEDS; i++)
	{
		seedx[i] = (int)((double)rand() / (double)RAND_MAX*(double)(sx - 1));
		seedy[i] = (int)((double)rand() / (double)RAND_MAX*(double)(sy - 1));
	}
	int *nseedx = new int[NSEEDS], *nseedy = new int[NSEEDS];
	for (int k = 0; k < 500; k++)
	{
		int seedx1 = 0, seedy1 = 0, seedx2 = 0, seedy2 = 0, seedx3 = 0, seedy3 = 0;
		for (int i = 0; i < NSEEDS; i++)
		{
			double distmin1 = 100000.0; int ind1 = 0;
			double distmin2 = 100000.0; int ind2 = 0;
			double distmin3 = 100000.0; int ind3 = 0;
			for (int j = 0; j < NSEEDS; j++)
			{
				if (i != j)
				{
					for (int ii = -1; ii <= 1; ii++) for (int jj = -1; jj <= 1; jj++)
					{
						double dist = sqrt((double)((seedx[i] - seedx[j] - ii*(sx))*(seedx[i] - seedx[j] - ii*(sx)) + (seedy[i] - seedy[j] - jj*(sy))*(seedy[i] - seedy[j] - jj*(sy))));
						if (dist < distmin1) { distmin3 = distmin2; ind3 = ind2; seedx3 = seedx2; seedy3 = seedy2; distmin2 = distmin1; ind2 = ind1; seedx2 = seedx1; seedy2 = seedy1; distmin1 = dist; ind1 = j; seedx1 = seedx[j] + ii*(sx); seedy1 = seedy[j] + jj*(sy); }
						else if (dist < distmin2) { distmin3 = distmin2; ind3 = ind2; seedx3 = seedx2; seedy3 = seedy2; distmin2 = dist; ind2 = j; seedx2 = seedx[j] + ii*(sx); seedy2 = seedy[j] + jj*(sy); }
						else if (dist < distmin3) { distmin3 = dist; ind3 = j; seedx3 = seedx[j] + ii*(sx); seedy3 = seedy[j] + jj*(sy); }
					}
				}
			}
			int dx1 = seedx1 - seedx[i];
			int dy1 = seedy1 - seedy[i];
			double no1 = sqrt((double)(dx1*dx1 + dy1*dy1)); if (no1 == 0.0) no1 = 1.0;
			int dx2 = seedx2 - seedx[i];
			int dy2 = seedy2 - seedy[i];
			double no2 = sqrt((double)(dx2*dx2 + dy2*dy2)); if (no2 == 0.0) no2 = 1.0;
			int dx3 = seedx3 - seedx[i];
			int dy3 = seedy3 - seedy[i];
			double no3 = sqrt((double)(dx3*dx3 + dy3*dy3)); if (no3 == 0.0) no3 = 1.0;
			double dirx = (double)dx1 / no1 / no1 + (double)dx2 / no2 / no2 + (double)dx3 / no3 / no3;
			double diry = (double)dy1 / no1 / no1 + (double)dy2 / no2 / no2 + (double)dy3 / no3 / no3;
			double no = sqrt(dirx*dirx + diry*diry); if (no == 0.0) no = 1.0;
			nseedx[i] = seedx[i] - (int)(dirx / no + 0.5); if (nseedx[i] < 0) nseedx[i] += sx; else if (nseedx[i] >= sx) nseedx[i] -= sx;
			nseedy[i] = seedy[i] - (int)(diry / no + 0.5); if (nseedy[i] < 0) nseedy[i] += sy; else if (nseedy[i] >= sy) nseedy[i] -= sy;
		}
		for (int i = 0; i < NSEEDS; i++) { seedx[i] = nseedx[i]; seedy[i] = nseedy[i]; }
	}
	delete[] nseedx;
	delete[] nseedy;
}

#define WITH_BLEND
#define WANG
#define THICKNESS 4



static const double TAMPX = 40.0;
static const double TAMPY = 40.0;
static const int NPATCH = 14;
static const double COL_THRESH = 30.0; // =0.0 no blending (=100.0 strong blending) on patch borders
static const int PYRRES = 128;
static const int MAXLEVELS = 5;
static const int NCOLORS = 8;
static const int MAXFRAGSIZE = 10;
static const int MINFRAGSIZE = 16 * 16;
static const int MAX = 5000; // search MAX translations, in reduced size
static const int RADIUS = 10;
static const int MARGIN = 2;
static const int MAX_TP_LEVELS = 6;
static const int fragthresh = 50;



void biscalenoisepatchexg(const std::string& filename_source, const std::string& base_dir, int NCLUSTERS, int CONTENTS, int NSCALES, int NROTANGLES, float TSCALE, int SZ_MULT)
{
	std::string filename_plus_ext = ASTex::IO::remove_path(filename_source);
	std::string name_file = ASTex::IO::remove_ext(filename_plus_ext);


	int levels = 0, ratio = 1;
	int resx, resy;

	hvPictRGB<unsigned char> pioriginal;
	ASTex::load_to_hv(filename_source, pioriginal);

	hvProcWPhase proc(base_dir + name_file + "_filtered_S.png", 1, false, false, false);

	proc.doTurbulence();
	hvPictRGB<unsigned char> turb(pioriginal.sizeX(), pioriginal.sizeY(), hvColRGB<unsigned char>(0));
	for (int i = 0; i < turb.sizeX(); i++) for (int j = 0; j < turb.sizeY(); j++)
	{
		double ppx = (double)i, ppy = (double)j;
		double tt = proc.turbulence(TSCALE*ppx, TSCALE*ppy);
		// Perlin turbulence
		if (tt < -1.0) tt = -1.0; else if (tt > 1.0) tt = 1.0;
		turb.update(i, j, hvColRGB<unsigned char>((unsigned char)(255.0 / 2.0*(tt + 1.0))));
	}

	// extract tile from wang

	hvPictRGB<unsigned char> wang, wangp, wangorig;

	ASTex::load_to_hv(base_dir + name_file + "_wang.png", wangorig);
	ASTex::load_to_hv(base_dir + name_file + "_filtw_wang.png", wang);
	ASTex::load_to_hv(base_dir + name_file + "_filt_wang.png", wangp);

	hvPictRGB<unsigned char> tile(wang.sizeX() / 2, wang.sizeY() / 2, hvColRGB<unsigned char>(0));
	hvPictRGB<unsigned char> tilep(wang.sizeX() / 2, wang.sizeY() / 2, hvColRGB<unsigned char>(0));
	hvPictRGB<unsigned char> tileorig(wang.sizeX() / 2, wang.sizeY() / 2, hvColRGB<unsigned char>(0));

	for (int i = 0; i < tile.sizeX() / 2; i++) for (int j = 0; j < tile.sizeY() / 2; j++)
	{
		tile.update(i, j, wang.get(i + tile.sizeX() / 2, j + tile.sizeY() / 2));
		tile.update(i + tile.sizeX() / 2, j, wang.get(i + 3 * tile.sizeX() / 2, j + tile.sizeY() / 2));
		tile.update(i, j + tile.sizeY() / 2, wang.get(i + tile.sizeX() / 2, j + 3 * tile.sizeY() / 2));
		tile.update(i + tile.sizeX() / 2, j + tile.sizeY() / 2, wang.get(i + 3 * tile.sizeX() / 2, j + 3 * tile.sizeY() / 2));

		tilep.update(i, j, wangp.get(i + tile.sizeX() / 2, j + tile.sizeY() / 2));
		tilep.update(i + tile.sizeX() / 2, j, wangp.get(i + 3 * tile.sizeX() / 2, j + tile.sizeY() / 2));
		tilep.update(i, j + tile.sizeY() / 2, wangp.get(i + tile.sizeX() / 2, j + 3 * tile.sizeY() / 2));
		tilep.update(i + tile.sizeX() / 2, j + tile.sizeY() / 2, wangp.get(i + 3 * tile.sizeX() / 2, j + 3 * tile.sizeY() / 2));

		tileorig.update(i, j, wangorig.get(i + tile.sizeX() / 2, j + tile.sizeY() / 2));
		tileorig.update(i + tile.sizeX() / 2, j, wangorig.get(i + 3 * tile.sizeX() / 2, j + tile.sizeY() / 2));
		tileorig.update(i, j + tile.sizeY() / 2, wangorig.get(i + tile.sizeX() / 2, j + 3 * tile.sizeY() / 2));
		tileorig.update(i + tile.sizeX() / 2, j + tile.sizeY() / 2, wangorig.get(i + 3 * tile.sizeX() / 2, j + 3 * tile.sizeY() / 2));
	}

	std::vector<hvPictRGB<unsigned char>> noisei(MAX_TP_LEVELS);
	for (int ind = 0; ind < NCLUSTERS; ind++)
	{
		ASTex::load_to_hv(base_dir + name_file + "_noise_synth" + std::to_string(ind) + ".png", noisei[ind]);
	}
	std::vector<hvPictRGB<unsigned char>> wangnoiseimask(MAX_TP_LEVELS);
	for (int ind = 0; ind < NCLUSTERS; ind++)
	{
		ASTex::load_to_hv(base_dir + name_file + "_noise_cl" + std::to_string(ind) + "_wang.png", wangnoiseimask[ind]);
	}
	std::vector< hvPictRGB<unsigned char>> noisetile(MAX_TP_LEVELS);
	std::vector<hvPict<unsigned char>> noisetileorig(MAX_TP_LEVELS);
	for (int ind = 0; ind < NCLUSTERS; ind++)
	{
		noisetile[ind].reset(wang.sizeX() / 2, wang.sizeY() / 2, hvColRGB<unsigned char>(0));
		for (int i = 0; i < tile.sizeX() / 2; i++) for (int j = 0; j < tile.sizeY() / 2; j++)
		{
			noisetile[ind].update(i, j, wangnoiseimask[ind].get(i + tile.sizeX() / 2, j + tile.sizeY() / 2));
			noisetile[ind].update(i + tile.sizeX() / 2, j, wangnoiseimask[ind].get(i + 3 * tile.sizeX() / 2, j + tile.sizeY() / 2));
			noisetile[ind].update(i, j + tile.sizeY() / 2, wangnoiseimask[ind].get(i + tile.sizeX() / 2, j + 3 * tile.sizeY() / 2));
			noisetile[ind].update(i + tile.sizeX() / 2, j + tile.sizeY() / 2, wangnoiseimask[ind].get(i + 3 * tile.sizeX() / 2, j + 3 * tile.sizeY() / 2));
		}
		hvBitmap noisem(wang.sizeX() / 2, wang.sizeY() / 2, false);
		for (int i = 0; i < tile.sizeX(); i++) for (int j = 0; j < tile.sizeY(); j++)
		{
			if (noisetile[ind].get(i, j).RED() > 0) noisem.set(i, j, true);
		}
		noisem.dilatation(3, 3); noisem.dilatation(3, 3); noisem.dilatation(3, 3);
		noisem.dilatation(3, 3); noisem.dilatation(3, 3); noisem.dilatation(3, 3);
		hvPict<unsigned char> noisempi(noisem, 4, 255);
		noisetileorig[ind].reset(wang.sizeX() / 2, wang.sizeY() / 2, 0);
		for (int i = 0; i < tile.sizeX(); i++) for (int j = 0; j < tile.sizeY(); j++)
		{
			noisetileorig[ind].update(i, j, noisempi.get(i, j));
		}
	}
	hvPictRGB<unsigned char> ptile[MAXLEVELS];
	if (tile.sizeX() > PYRRES || tile.sizeY() > PYRRES)
		do {
			ptile[levels].shrink((levels == 0 ? &tile : ptile + (levels - 1)));
			levels++;
		} while (levels < MAXLEVELS && (ptile[levels - 1].sizeX() > PYRRES || ptile[levels - 1].sizeY() > PYRRES));
		ratio = 1 << levels;
		resx = (levels > 0 ? ptile[levels - 1].sizeX() : tile.sizeX());
		resy = (levels > 0 ? ptile[levels - 1].sizeY() : tile.sizeY());

		hvPictRGB<unsigned char> weight(tile.sizeX(), tile.sizeY(), hvColRGB<unsigned char>(0));
		for (int i = 0; i < tile.sizeX(); i++) for (int j = 0; j < tile.sizeY(); j++)
		{
			double rr = 0.0, gg = 0.0, bb = 0.0, sum = 0.0;
			for (int ii = -THICKNESS; ii <= THICKNESS; ii++) for (int jj = -THICKNESS; jj <= THICKNESS; jj++)
			{
				double ff = exp(-(ii*ii + jj*jj) / 8.0);
				sum += ff;
				int xx = i + ii; if (xx < 0) xx += tile.sizeX(); else if (xx >= tile.sizeX()) xx -= tile.sizeX();
				int yy = j + jj; if (yy < 0) yy += tile.sizeY(); else if (yy >= tile.sizeY()) yy -= tile.sizeY();
				hvColRGB<unsigned char> cc; cc = tilep.get(xx, yy);
				rr += ff*(double)cc.RED();
				gg += ff*(double)cc.GREEN();
				bb += ff*(double)cc.BLUE();
			}
			hvColRGB<unsigned char> colg((unsigned char)(rr / sum), (unsigned char)(gg / sum), (unsigned char)(bb / sum));
			hvColRGB<unsigned char> col = tilep.get(i, j);
			int vv = (int)sqrt(colg.squaredDifference(col)) * 4; if (vv > 255) vv = 255;
			weight.update(i, j, hvColRGB<unsigned char>((unsigned char)vv));
		}

		Fragmentation::Fragmenter f(&tile, true);
		f.fragment(tile.sizeX()*tile.sizeY() / MAXFRAGSIZE / MAXFRAGSIZE, fragthresh);
		f.cleanup(MINFRAGSIZE);


		hvQPictRGB<unsigned char, 64> wangq(tile, levels > 0 ? ptile[levels - 1] : tile, NCOLORS);
		hvPictRGB<unsigned char> lob;
		wangq.convert(lob, 1);

		hvPictRGB<unsigned char> tilesal(tile.sizeX(), tile.sizeY(), hvColRGB<unsigned char>(0));
		for (int i = 0; i < tilesal.sizeX(); i++) for (int j = 0; j < tilesal.sizeY(); j++)
		{
			int ind = wangq.closest(tilep.get(i, j));
			tilesal.update(i, j, wangq.getColor(ind));
		}


		// make the fragment map "fragpf" : X=frag index, Y=delta_x, Z=delta_y and W=corner (0,1,2 or 3)
		////////////////////////////////////////////////////
		hvPict<hvVec4<int> > fragpf(tile.sizeX(), tile.sizeY(), hvVec4<int>(0));
		hvBitmap yet(tile.sizeX(), tile.sizeY(), false);
		hvBitmap maskl(tile.sizeX(), tile.sizeY(), false);
		hvBitmap maskr(tile.sizeX(), tile.sizeY(), false);
		hvBitmap maskt(tile.sizeX(), tile.sizeY(), false);
		hvBitmap maskb(tile.sizeX(), tile.sizeY(), false);
		for (int i = 0; i < f.nFragments(); i++)
		{
			bool bl = false, br = false, bt = false, bb = false;
			hvVec2<int> min, max;
			int meanx = 0, meany = 0;
			int count = f.fragmentSize(i);
			for (int j = 0; j < count; j++)
			{
				int x, y; f.pixelFragment(i, j, x, y);
				if (j == 0) { min = hvVec2<int>(x, y); max == hvVec2<int>(x, y); }
				else
				{
					min.keepMin(min, hvVec2<int>(x, y));
					max.keepMax(max, hvVec2<int>(x, y));
				}
				if (x == 0) bl = true;
				if (x == tile.sizeX() - 1) br = true;
				if (y == 0) bb = true;
				if (y == tile.sizeY() - 1) bt = true;
				meanx += x; meany += y;
			}
			meanx /= count; meany /= count;
			if ((bl&&br) || (bt&&bb))
			{
				int ii, jj, k;
				yet.clear(true);
				for (int j = 0; j < count; j++)
				{
					int x, y; f.pixelFragment(i, j, x, y);
					yet.set(x, y, false);
				}

				meanx = 0; meany = 0;
				if (bl && br && !(bt&&bb)) // frag has only vertical border
				{
					hvVec2<int> bminl, bmaxl;
					hvVec2<int> bminr, bmaxr;
					maskl.clear(false);
					maskr.clear(false);
					bool first = true;
					int cc = 0;
					for (k = 0; k < tile.sizeY(); k++) if (!yet.get(0, k)) { bminl = hvVec2<int>(0, k); bmaxl = hvVec2<int>(0, k); }
					for (k = 0; k < tile.sizeY(); k++) if (!yet.get(0, k)) { yet.seedfill(0, k, maskl, bminl, bmaxl); }
					if (bminl.X() == 0 && bmaxl.X() == tile.sizeX() - 1) hvFatal("fragment over opposite vertical borders");
					//if (k==tile.sizeY()) hvFatal("could not start seed (a)");
					for (int ii = bminl.X(); ii <= bmaxl.X(); ii++) for (int jj = bminl.Y(); jj <= bmaxl.Y(); jj++)
					{
						if (maskl.get(ii, jj)) { cc++; meanx += (ii + tile.sizeX()); meany += jj; }
					}
					for (k = 0; k < tile.sizeY(); k++) if (!yet.get(tile.sizeX() - 1, k)) { bminr = hvVec2<int>(tile.sizeX() - 1, k); bmaxr = hvVec2<int>(tile.sizeX() - 1, k); }
					for (k = 0; k < tile.sizeY(); k++) if (!yet.get(tile.sizeX() - 1, k)) { yet.seedfill(tile.sizeX() - 1, k, maskr, bminr, bmaxr); }
					if (bminr.X() == 0 && bmaxr.X() == tile.sizeX() - 1) hvFatal("fragment over opposite vertical borders");
					//if (k==tile.sizeY()) hvFatal("could not start seed (b)");
					for (int ii = bminr.X(); ii <= bmaxr.X(); ii++) for (int jj = bminr.Y(); jj <= bmaxr.Y(); jj++)
					{
						if (maskr.get(ii, jj)) { cc++; meanx += ii; meany += jj; }
					}
					for (int ii = 0; ii < tile.sizeX(); ii++) for (int jj = 0; jj < tile.sizeY(); jj++) if (!yet.get(ii, jj)) hvFatal("fragment has more than one connected component");
					if (cc != count) hvFatal("fragment has more than one connected component");
					meanx /= count; meany /= count;
					for (int j = 0; j < count; j++)
					{
						int x, y; f.pixelFragment(i, j, x, y);
						fragpf.update(x, y, hvVec4<int>(i, meanx, meany, maskl.get(x, y) ? 1 : 0));
					}
					//gets(buff);
				}
				else if (bt && bb && !(bl&&br)) // frag has only horizontal border
				{
					hvVec2<int> bminb, bmaxb;
					hvVec2<int> bmint, bmaxt;
					maskb.clear(false);
					maskt.clear(false);
					bool first = true;
					int cc = 0;
					for (k = 0; k < tile.sizeX(); k++) if (!yet.get(k, 0)) { bminb = hvVec2<int>(k, 0); bmaxb = hvVec2<int>(k, 0); break; }
					for (k = 0; k < tile.sizeX(); k++) if (!yet.get(k, 0)) { yet.seedfill(k, 0, maskb, bminb, bmaxb); }
					if (bminb.Y() == 0 && bmaxb.Y() == tile.sizeY() - 1) hvFatal("fragment over opposite horizontal borders");

					for (int ii = bminb.X(); ii <= bmaxb.X(); ii++) for (int jj = bminb.Y(); jj <= bmaxb.Y(); jj++)
					{
						if (maskb.get(ii, jj)) { cc++; meanx += ii; meany += jj + tile.sizeY(); }
					}
					for (k = 0; k < tile.sizeX(); k++) if (!yet.get(k, tile.sizeY() - 1)) { bmint = hvVec2<int>(k, tile.sizeY() - 1); bmaxt = hvVec2<int>(k, tile.sizeY() - 1); break; }
					for (k = 0; k < tile.sizeX(); k++) if (!yet.get(k, tile.sizeY() - 1)) { yet.seedfill(k, tile.sizeY() - 1, maskt, bmint, bmaxt); }
					if (bmint.Y() == 0 && bmaxt.Y() == tile.sizeY() - 1) hvFatal("fragment over opposite horizontal borders");

					for (int ii = bmint.X(); ii <= bmaxt.X(); ii++) for (int jj = bmint.Y(); jj <= bmaxt.Y(); jj++)
					{
						if (maskt.get(ii, jj)) { cc++; meanx += ii; meany += jj; }
					}

					for (int ii = 0; ii < tile.sizeX(); ii++) for (int jj = 0; jj < tile.sizeY(); jj++) if (!yet.get(ii, jj)) hvFatal("fragment has more than one connected component");
					if (cc != count) hvFatal("fragment has more than one connected component");
					meanx /= count; meany /= count;
					for (int j = 0; j < count; j++)
					{
						int x, y; f.pixelFragment(i, j, x, y);
						fragpf.update(x, y, hvVec4<int>(i, meanx, meany, maskb.get(x, y) ? 2 : 0));
					}
					//gets(buff);
				}
				else
				{
					maskl.clear(false);
					maskr.clear(false);
					maskb.clear(false);
					maskt.clear(false);
					bool finish = true;
					for (k = 0; k < tile.sizeX() / 2 - 1 && finish; k++)
					{
						int xx;
						for (xx = 0; xx <= k && finish; xx++) if (!yet.get(xx, k)) { yet.seedfillCorners(xx, k, maskl, maskr, maskb, maskt, 3, meanx, meany); finish = false; }
						for (xx = 0; xx <= k && finish; xx++) if (!yet.get(k, xx)) { yet.seedfillCorners(k, xx, maskl, maskr, maskb, maskt, 3, meanx, meany); finish = false; }
						for (xx = 0; xx <= k && finish; xx++) if (!yet.get(tile.sizeX() - 1 - xx, k)) { yet.seedfillCorners(tile.sizeX() - 1 - xx, k, maskl, maskr, maskb, maskt, 2, meanx, meany); finish = false; }
						for (xx = 0; xx <= k && finish; xx++) if (!yet.get(tile.sizeX() - 1 - k, xx)) { yet.seedfillCorners(tile.sizeX() - 1 - k, xx, maskl, maskr, maskb, maskt, 2, meanx, meany); finish = false; }
						for (xx = 0; xx <= k && finish; xx++) if (!yet.get(xx, tile.sizeY() - 1 - k)) { yet.seedfillCorners(xx, tile.sizeY() - 1 - k, maskl, maskr, maskb, maskt, 1, meanx, meany); finish = false; }
						for (xx = 0; xx <= k && finish; xx++) if (!yet.get(k, tile.sizeY() - 1 - xx)) { yet.seedfillCorners(k, tile.sizeY() - 1 - xx, maskl, maskr, maskb, maskt, 1, meanx, meany); finish = false; }
						for (xx = 0; xx <= k && finish; xx++) if (!yet.get(tile.sizeX() - 1 - xx, tile.sizeY() - 1 - k)) { yet.seedfillCorners(tile.sizeX() - 1 - xx, tile.sizeY() - 1 - k, maskl, maskr, maskb, maskt, 0, meanx, meany); finish = false; }
						for (xx = 0; xx <= k && finish; xx++) if (!yet.get(tile.sizeX() - 1 - k, tile.sizeY() - 1 - xx)) { yet.seedfillCorners(tile.sizeX() - 1 - k, tile.sizeY() - 1 - xx, maskl, maskr, maskb, maskt, 0, meanx, meany); finish = false; }
					}
					for (int ii = 0; ii < tile.sizeX(); ii++) for (int jj = 0; jj < tile.sizeY(); jj++) if (!yet.get(ii, jj)) hvFatal("fragment has more than one connected component");
					meanx /= count; meany /= count;
					for (int j = 0; j < count; j++)
					{
						int x, y; f.pixelFragment(i, j, x, y);
						int lab = 0;
						if (maskl.get(x, y)) lab = 3;
						else if (maskr.get(x, y)) lab = 2;
						else if (maskb.get(x, y)) lab = 1;
						fragpf.update(x, y, hvVec4<int>(i, meanx, meany, lab));
					}
					//gets(buff);
				}
			}
			else for (int j = 0; j < count; j++)
			{
				int x, y; f.pixelFragment(i, j, x, y);
				fragpf.update(x, y, hvVec4<int>(i, meanx, meany, 0));
			}
		}

		hvPictRGB<unsigned char> pifrag(tile.sizeX(), tile.sizeY(), hvColRGB<unsigned char>(0));
		for (int i = 0; i < tile.sizeX(); i++)
		{
			for (int j = 0; j < tile.sizeY(); j++)
			{
				hvVec4<int> col = fragpf.get(i, j);
				unsigned char px = (unsigned char)((double)((double)col.Y() - (double)(col.W() == 3 || col.W() == 1 ? tile.sizeX() : 0) - (double)i) / 512.0*127.0 + 128.0);
				unsigned char py = (unsigned char)((double)((double)col.Z() - (double)(col.W() == 3 || col.W() == 2 ? tile.sizeY() : 0) - (double)j) / 512.0*127.0 + 128.0);
				pifrag.update(i, j, hvColRGB<unsigned char>(px, py, 128));
			}
		}

		// create patches: a patch = set of fragments
		//////////////////////////////////////////////
		int ii, jj;
		hvPict<int> frags(tile.sizeX(), tile.sizeY(), 0);
		for (int i = 0; i < f.nFragments(); i++)
		{
			int count = f.fragmentSize(i);
			for (int j = 0; j < count; j++)
			{
				int x, y; f.pixelFragment(i, j, x, y);
				frags.update(x, y, i);
			}
		}
		int count = f.nFragments();
		hvSortedList<int> **fneighbors = new hvSortedList<int>*[count];
		bool *stillfree = new bool[count];
		for (int i = 0; i < count; i++) { fneighbors[i] = new hvSortedList<int>(500); stillfree[i] = true; }
		for (int i = 0; i < f.nFragments(); i++)
		{
			int count = f.fragmentSize(i);
			for (int j = 0; j < count; j++)
			{
				int x, y; f.pixelFragment(i, j, x, y);
				int k;
				int nnx[4] = { 0,0,-1,1 }, nny[4] = { -1,1,0,0 };
				for (int k = 0; k < 4; k++)
				{
					ii = nnx[k]; jj = nny[k];
					int xx = x + ii; if (xx < 0) xx += tile.sizeX(); else if (xx >= tile.sizeX()) xx -= tile.sizeX();
					int yy = y + jj; if (yy < 0) yy += tile.sizeY(); else if (yy >= tile.sizeY()) yy -= tile.sizeY();
					int ind = frags.get(xx, yy);
					if (ind != i) fneighbors[i]->pushSortedSet(ind);
				}
			}
		}


		int seedx[NPATCH], seedy[NPATCH];
		Poisson(tile.sizeX(), tile.sizeY(), NPATCH, seedx, seedy);
		hvSortedList<int> *patchfrags[NPATCH], *candidates[NPATCH];
		for (int i = 0; i < NPATCH; i++)
		{
			patchfrags[i] = new hvSortedList<int>(count); patchfrags[i]->clear();
			candidates[i] = new hvSortedList<int>(count); candidates[i]->clear();
			int ind = frags.get(seedx[i], seedy[i]);
			if (!stillfree[ind])
			{
				bool found = true;
				int aa = ind;
				for (int j = 0; j < fneighbors[ind]->size() && found; j++) if (stillfree[fneighbors[ind]->at(j)]) { aa = fneighbors[ind]->at(j); found = false; }
				if (found) hvFatal("cannot find initial fragment for patch");
				ind = aa;
			}
			patchfrags[i]->push_back(ind);
			stillfree[ind] = false;
			for (int j = 0; j < fneighbors[ind]->size(); j++) if (stillfree[fneighbors[ind]->at(j)]) candidates[i]->push_back(fneighbors[ind]->at(j));
		}
		bool remains = true;
		do {
			remains = false;
			for (int i = 0; i < NPATCH; i++)
			{
				if (candidates[i]->size() > 0)
				{
					int fragind = candidates[i]->at(0);
					candidates[i]->erase(candidates[i]->begin());
					patchfrags[i]->push_back(fragind);
					stillfree[fragind] = false;
					remains = true;
					for (int j = 0; j < fneighbors[fragind]->size(); j++)
					{
						int nind = fneighbors[fragind]->at(j);
						if (stillfree[nind] && candidates[i]->search(nind) == -1) candidates[i]->push_back(nind);
					}
					for (int j = 0; j < NPATCH; j++) if (j != i)
					{
						bool removed = true;
						for (int k = 0; k < candidates[j]->size() && removed; k++) if (candidates[j]->at(k) == fragind) { removed = false; candidates[j]->erase(candidates[j]->begin() + k); }
					}
				}
			}
		} while (remains);


		hvPictRGB<unsigned char> patchmaprgb(tile.sizeX(), tile.sizeY(), hvColRGB<unsigned char>());
		hvPict<int> patchmap(tile.sizeX(), tile.sizeY(), 0);
		for (int i = 0; i < NPATCH; i++)
		{
			for (int j = 0; j < patchfrags[i]->size(); j++)
			{
				int fragind = patchfrags[i]->at(j);
				for (int ii = 0; ii < f.fragmentSize(fragind); ii++)
				{
					int x, y; f.pixelFragment(fragind, ii, x, y);
					patchmaprgb.update(x, y, hvColRGB<unsigned char>(((i + 11) * 33) % 255, ((i + 54) * 11) % 255, ((i + 33) * 67) % 255));
					patchmap.update(x, y, i);
				}
			}
		}

		for (int i = 0; i < count; i++) { delete fneighbors[i]; }
		delete[] fneighbors;
		delete[] stillfree;


		// create patchmap "patchpf" like fragment map
		/////////////////////////////////////////////
		std::vector<hvVec2<int> > *patchpixels[NPATCH];
		hvPict<hvVec4<int> > patchpf(tile.sizeX(), tile.sizeY(), hvVec4<int>(0));
		yet.clear(false); maskl.clear(false); maskr.clear(false); maskt.clear(false); maskb.clear(false);
		for (int i = 0; i < NPATCH; i++)
		{
			bool bl = false, br = false, bt = false, bb = false;
			hvVec2<int> min, max;
			int meanx = 0, meany = 0;
			int count = 0;
			for (int j = 0; j < patchfrags[i]->size(); j++)
			{
				int k;
				int fragind = patchfrags[i]->at(j);
				count += f.fragmentSize(fragind);
			}
			patchpixels[i] = new std::vector<hvVec2<int> >(count);
			yet.clear(true);
			for (int j = 0; j < patchfrags[i]->size(); j++)
			{
				int k;
				int fragind = patchfrags[i]->at(j);
				for (k = 0; k < f.fragmentSize(fragind); k++)
				{
					int x, y; f.pixelFragment(fragind, k, x, y);
					yet.set(x, y, false);
					if (j == 0 && k == 0) { min = hvVec2<int>(x, y); max == hvVec2<int>(x, y); }
					else
					{
						min.keepMin(min, hvVec2<int>(x, y));
						max.keepMax(max, hvVec2<int>(x, y));
					}
					if (x == 0) bl = true;
					if (x == tile.sizeX() - 1) br = true;
					if (y == 0) bb = true;
					if (y == tile.sizeY() - 1) bt = true;
					meanx += x; meany += y;
				}
			}
			hvPictRGB<unsigned char> patchcurr(yet, hvColRGB<unsigned char>(0), hvColRGB<unsigned char>(255));

			meanx /= count; meany /= count;
			if ((bl&&br) || (bt&&bb)) // patch is on a border
			{
				int ii, jj, k;
				yet.clear(true);
				for (int j = 0; j < patchfrags[i]->size(); j++)
				{
					int fragind = patchfrags[i]->at(j);
					for (k = 0; k < f.fragmentSize(fragind); k++)
					{
						int x, y; f.pixelFragment(fragind, k, x, y);
						yet.set(x, y, false);
					}
				}
				meanx = 0; meany = 0;
				if (bl && br && !(bt&&bb)) // frag has only vertical border
				{
					hvVec2<int> bminl, bmaxl;
					hvVec2<int> bminr, bmaxr;
					maskl.clear(false);
					maskr.clear(false);
					bool first = true;
					int cc = 0;
					for (k = 0; k < tile.sizeY(); k++) if (!yet.get(0, k)) { bminl = hvVec2<int>(0, k); bmaxl = hvVec2<int>(0, k); }
					for (k = 0; k < tile.sizeY(); k++) if (!yet.get(0, k)) { yet.seedfill(0, k, maskl, bminl, bmaxl); }
					if (bminl.X() == 0 && bmaxl.X() == tile.sizeX() - 1) hvFatal("patch over opposite vertical borders");

					for (int ii = bminl.X(); ii <= bmaxl.X(); ii++) for (int jj = bminl.Y(); jj <= bmaxl.Y(); jj++)
					{
						if (maskl.get(ii, jj)) { cc++; meanx += (ii + tile.sizeX()); meany += jj; }
					}
					for (k = 0; k < tile.sizeY(); k++) if (!yet.get(tile.sizeX() - 1, k)) { bminr = hvVec2<int>(tile.sizeX() - 1, k); bmaxr = hvVec2<int>(tile.sizeX() - 1, k); }
					for (k = 0; k < tile.sizeY(); k++) if (!yet.get(tile.sizeX() - 1, k)) { yet.seedfill(tile.sizeX() - 1, k, maskr, bminr, bmaxr); }
					if (bminr.X() == 0 && bmaxr.X() == tile.sizeX() - 1) hvFatal("patch over opposite vertical borders");

					for (int ii = bminr.X(); ii <= bmaxr.X(); ii++) for (int jj = bminr.Y(); jj <= bmaxr.Y(); jj++)
					{
						if (maskr.get(ii, jj)) { cc++; meanx += ii; meany += jj; }
					}

					meanx /= count; meany /= count;
					for (int j = 0; j < patchfrags[i]->size(); j++)
					{
						int fragind = patchfrags[i]->at(j);
						for (k = 0; k < f.fragmentSize(fragind); k++)
						{
							int x, y; f.pixelFragment(fragind, k, x, y);
							patchpf.update(x, y, hvVec4<int>(i, meanx, meany, maskl.get(x, y) ? 1 : 0));
							patchpixels[i]->push_back(hvVec2<int>(x + maskl.get(x, y) ? tile.sizeX() : 0, y));
						}
					}
					//gets(buff);
				}
				else if (bt && bb && !(bl&&br)) // frag has only horizontal border
				{
					hvVec2<int> bminb, bmaxb;
					hvVec2<int> bmint, bmaxt;
					maskb.clear(false);
					maskt.clear(false);
					bool first = true;
					int cc = 0;
					for (k = 0; k < tile.sizeX(); k++) if (!yet.get(k, 0)) { bminb = hvVec2<int>(k, 0); bmaxb = hvVec2<int>(k, 0); break; }
					for (k = 0; k < tile.sizeX(); k++) if (!yet.get(k, 0)) { yet.seedfill(k, 0, maskb, bminb, bmaxb); }
					if (bminb.Y() == 0 && bmaxb.Y() == tile.sizeY() - 1) hvFatal("patch over opposite horizontal borders");
					for (int ii = bminb.X(); ii <= bmaxb.X(); ii++) for (int jj = bminb.Y(); jj <= bmaxb.Y(); jj++)
					{
						if (maskb.get(ii, jj)) { cc++; meanx += ii; meany += jj + tile.sizeY(); }
					}
					for (k = 0; k < tile.sizeX(); k++) if (!yet.get(k, tile.sizeY() - 1)) { bmint = hvVec2<int>(k, tile.sizeY() - 1); bmaxt = hvVec2<int>(k, tile.sizeY() - 1); break; }
					for (k = 0; k < tile.sizeX(); k++) if (!yet.get(k, tile.sizeY() - 1)) { yet.seedfill(k, tile.sizeY() - 1, maskt, bmint, bmaxt); }
					if (bmint.Y() == 0 && bmaxt.Y() == tile.sizeY() - 1) hvFatal("patch over opposite horizontal borders");

					for (int ii = bmint.X(); ii <= bmaxt.X(); ii++) for (int jj = bmint.Y(); jj <= bmaxt.Y(); jj++)
					{
						if (maskt.get(ii, jj)) { cc++; meanx += ii; meany += jj; }
					}

					meanx /= count; meany /= count;
					for (int j = 0; j < patchfrags[i]->size(); j++)
					{
						int fragind = patchfrags[i]->at(j);
						for (k = 0; k < f.fragmentSize(fragind); k++)
						{
							int x, y; f.pixelFragment(fragind, k, x, y);
							patchpf.update(x, y, hvVec4<int>(i, meanx, meany, maskb.get(x, y) ? 2 : 0));
							patchpixels[i]->push_back(hvVec2<int>(x, y + maskb.get(x, y) ? tile.sizeY() : 0));
						}
					}
				}
				else
				{
					maskl.clear(false);
					maskr.clear(false);
					maskb.clear(false);
					maskt.clear(false);
					bool finish = true;
					for (k = 0; k < tile.sizeX() / 2 - 1 && finish; k++)
					{
						int xx;
						for (xx = 0; xx <= k && finish; xx++) if (!yet.get(xx, k)) { yet.seedfillCorners(xx, k, maskl, maskr, maskb, maskt, 3, meanx, meany); finish = false; }
						for (xx = 0; xx <= k && finish; xx++) if (!yet.get(k, xx)) { yet.seedfillCorners(k, xx, maskl, maskr, maskb, maskt, 3, meanx, meany); finish = false; }
						for (xx = 0; xx <= k && finish; xx++) if (!yet.get(tile.sizeX() - 1 - xx, k)) { yet.seedfillCorners(tile.sizeX() - 1 - xx, k, maskl, maskr, maskb, maskt, 2, meanx, meany); finish = false; }
						for (xx = 0; xx <= k && finish; xx++) if (!yet.get(tile.sizeX() - 1 - k, xx)) { yet.seedfillCorners(tile.sizeX() - 1 - k, xx, maskl, maskr, maskb, maskt, 2, meanx, meany); finish = false; }
						for (xx = 0; xx <= k && finish; xx++) if (!yet.get(xx, tile.sizeY() - 1 - k)) { yet.seedfillCorners(xx, tile.sizeY() - 1 - k, maskl, maskr, maskb, maskt, 1, meanx, meany); finish = false; }
						for (xx = 0; xx <= k && finish; xx++) if (!yet.get(k, tile.sizeY() - 1 - xx)) { yet.seedfillCorners(k, tile.sizeY() - 1 - xx, maskl, maskr, maskb, maskt, 1, meanx, meany); finish = false; }
						for (xx = 0; xx <= k && finish; xx++) if (!yet.get(tile.sizeX() - 1 - xx, tile.sizeY() - 1 - k)) { yet.seedfillCorners(tile.sizeX() - 1 - xx, tile.sizeY() - 1 - k, maskl, maskr, maskb, maskt, 0, meanx, meany); finish = false; }
						for (xx = 0; xx <= k && finish; xx++) if (!yet.get(tile.sizeX() - 1 - k, tile.sizeY() - 1 - xx)) { yet.seedfillCorners(tile.sizeX() - 1 - k, tile.sizeY() - 1 - xx, maskl, maskr, maskb, maskt, 0, meanx, meany); finish = false; }
					}
					meanx /= count; meany /= count;
					for (int j = 0; j < patchfrags[i]->size(); j++)
					{
						int fragind = patchfrags[i]->at(j);
						for (k = 0; k < f.fragmentSize(fragind); k++)
						{
							int x, y; f.pixelFragment(fragind, k, x, y);
							int lab = 0;
							if (maskl.get(x, y)) lab = 3;
							else if (maskr.get(x, y)) lab = 2;
							else if (maskb.get(x, y)) lab = 1;
							patchpf.update(x, y, hvVec4<int>(i, meanx, meany, lab));
							patchpixels[i]->push_back(hvVec2<int>(x + ((lab == 3 || lab == 1) ? tile.sizeX() : 0), y + ((lab == 3 || lab == 2) ? tile.sizeY() : 0)));
						}
					}
				}
			}
			else for (int j = 0; j < patchfrags[i]->size(); j++)
			{
				int k;
				int fragind = patchfrags[i]->at(j);
				for (k = 0; k < f.fragmentSize(fragind); k++)
				{
					int x, y; f.pixelFragment(fragind, k, x, y);
					patchpf.update(x, y, hvVec4<int>(i, meanx, meany, 0));
					patchpixels[i]->push_back(hvVec2<int>(x, y));
				}
			}
		}

		hvPictRGB<unsigned char> pipatch(tile.sizeX(), tile.sizeY(), hvColRGB<unsigned char>(0));
		for (int i = 0; i < tile.sizeX(); i++)
		{
			for (int j = 0; j < tile.sizeY(); j++)
			{
				hvVec4<int> col = patchpf.get(i, j);
				unsigned char px = (unsigned char)((double)((double)col.Y() - (double)(col.W() == 3 || col.W() == 1 ? tile.sizeX() : 0) - (double)i) / (double)tile.sizeX()*127.0 + 128.0);
				unsigned char py = (unsigned char)((double)((double)col.Z() - (double)(col.W() == 3 || col.W() == 2 ? tile.sizeY() : 0) - (double)j) / (double)tile.sizeY()*127.0 + 128.0);
				pipatch.update(i, j, hvColRGB<unsigned char>(px, py, 128));
			}
		}

		// make patch erosion map
		hvPictRGB<unsigned char> patcherod(tile.sizeX(), tile.sizeY(), 0);
		for (int i = 0; i < NPATCH; i++)
		{
			hvBitmap patchmask(tile.sizeX(), tile.sizeY(), false);
			for (int j = 0; j < patchfrags[i]->size(); j++)
			{
				int k;
				int fragind = patchfrags[i]->at(j);
				for (k = 0; k < f.fragmentSize(fragind); k++)
				{
					int u, v; f.pixelFragment(fragind, k, u, v);
					patchmask.set(u, v, true);
				}
			}
			hvPict<unsigned char> perod(true, patchmask, THICKNESS + 1, 255);
			for (int j = 0; j < patchfrags[i]->size(); j++)
			{
				int k;
				int fragind = patchfrags[i]->at(j);
				for (k = 0; k < f.fragmentSize(fragind); k++)
				{
					int u, v; f.pixelFragment(fragind, k, u, v);
					patcherod.update(u, v, hvColRGB<unsigned char>(perod.get(u, v)));
				}
			}
		}

		// find contents
		/////////////////////////////////////////////

		std::vector<std::vector<int>>    deltax(100, std::vector<int>(CONTENTS));
		std::vector<std::vector<int>>    deltay(100, std::vector<int>(CONTENTS));
		std::vector<std::vector<double>> scalex(100, std::vector<double>(CONTENTS));
		std::vector<std::vector<double>> scaley(100, std::vector<double>(CONTENTS));
		std::vector<std::vector<int>>    idscale(100, std::vector<int>(CONTENTS));
		std::vector<std::vector<int>>    idrot(100, std::vector<int>(CONTENTS));
		std::vector<std::vector<double>> errpatch(100, std::vector<double>(CONTENTS));
		int cx[100], cy[100];
		int ccx[100], ccy[100];

		double sfactorx[3] = { 1.0, 1.5, 2.0 };
		double sfactory[3] = { 1.0, 1.5, 2.0 };
		double rotangles[5] = { 0.0, M_PI / 2.0, M_PI, M_PI / 4.0,  3.0*M_PI / 4.0 };

		for (int npatch = 0; npatch < NPATCH; npatch++)
		{
			hvBitmap patchb(2 * resx, 2 * resy, false);
			hvBitmap feature(2 * resx, 2 * resy, false);
			for (int j = 0; j < patchfrags[npatch]->size(); j++)
			{
				int k;
				int fragind = patchfrags[npatch]->at(j);
				for (k = 0; k < f.fragmentSize(fragind); k++)
				{
					int u, v; f.pixelFragment(fragind, k, u, v);
					hvVec4<int> col = patchpf.get(u, v);
					int x = u / ratio + (col.W() == 3 || col.W() == 1 ? resx : 0);
					int y = v / ratio + (col.W() == 3 || col.W() == 2 ? resy : 0);
					patchb.set(x, y, true);
				}
			}
			patchb.fillholes(); feature = patchb;
			feature.erosion(3, 3); ~feature;
			patchb.dilatation(3, 3); patchb.dilatation(3, 3); //patch.fillholes();
																					   //~feature; feature.dilatation(3,3); feature.fillholes(); ~feature;
			patchb &= feature;

			std::vector<hvVec2<int> > border(10000); border.clear();
			int minx = 2 * resx, miny = 2 * resy, maxx = 0, maxy = 0;
			for (int i = 0; i < 2 * resx; i++)
			{
				for (int j = 0; j < 2 * resy; j++)
				{
					if (patchb.get(i, j))
					{
						border.push_back(hvVec2<int>(i, j));
						if (minx > i) minx = i;
						if (maxx < i) maxx = i;
						if (miny > j) miny = j;
						if (maxy < j) maxy = j;
					}
				}
			}

			cx[npatch] = (maxx + minx) / 2; cy[npatch] = (maxy + miny) / 2;

			std::vector<int> dminx(MAX);
			std::vector<int> dminy(MAX);
			std::vector<int> smin(MAX);
			std::vector<int> rotmin(MAX);
			std::vector<double> errormin(MAX, 1000000000.0);
			double maxerror = 0.0;

			hvPict<double> patcherr(resx, resy, 0.0);

			for (int kr = 0; kr < NROTANGLES; kr++)
				for (int ns = 0; ns < NSCALES; ns++)
					for (int j = 0; j < resy; j += 1)
#pragma omp parallel for
						for (int i = 0; i < resx; i += 1)
						{
							double errsum = 0.0;
							for (int k = 0; k < border.size(); k++)
							{
								hvVec2<int> pos = border.at(k);
								double theta = rotangles[kr];
								int bx = (int)((double)(pos.X() - cx[npatch])*sfactorx[ns] * cos(theta) - (double)(pos.Y() - cy[npatch])*sfactory[ns] * sin(theta)) + i;
								int by = (int)((double)(pos.X() - cx[npatch])*sfactorx[ns] * sin(theta) + (double)(pos.Y() - cy[npatch])*sfactory[ns] * cos(theta)) + j;
								while (bx < 0) bx += resx; while (bx >= resx) bx -= resx;
								while (by < 0) by += resy; while (by >= resy) by -= resy;
								hvColRGB<unsigned char> cola = (levels > 0 ? ptile[levels - 1].get(bx, by) : tile.get(bx, by));
								int ppx = pos.X(); ppx = ppx%resx;
								int ppy = pos.Y(); ppy = ppy%resy;
								hvColRGB<unsigned char> colb = (levels > 0 ? ptile[levels - 1].get(ppx, ppy) : tile.get(ppx, ppy));
								double mm = cola.squaredDifferenceNorm(255.0, colb);
								errsum += mm;
							}
							errsum /= (double)border.size();
							int k;
							for (k = 0; k < MAX; k++) if (errsum < errormin[k]) break;
							if (k < MAX)
							{
								for (int m = MAX - 2; m >= k; m--)
								{
									errormin[m + 1] = errormin[m];
									dminx[m + 1] = dminx[m]; dminy[m + 1] = dminy[m];
									smin[m + 1] = smin[m]; rotmin[m + 1] = rotmin[m];
								}
								errormin[k] = errsum; dminx[k] = i; dminy[k] = j; smin[k] = ns; rotmin[k] = kr;
							}
							patcherr.update(i, j, errsum);
							if (maxerror < errsum) maxerror = errsum;
						}

			hvPictRGB<unsigned char> patchexgerr(resx, resy, hvColRGB<unsigned char>());

			for (int j = 0; j < resy; j += 1)
#pragma omp parallel for
				for (int i = 0; i < resx; i += 1)
				{
					patchexgerr.update(i, j, hvColRGB<unsigned char>((unsigned char)(patcherr.get(i, j) / maxerror*255.0)));
				}

			deltax[npatch][0] = dminx[0]; deltay[npatch][0] = dminy[0];
			scalex[npatch][0] = sfactorx[smin[0]]; scaley[npatch][0] = sfactory[smin[0]];
			idscale[npatch][0] = smin[0];
			idrot[npatch][0] = rotmin[0];
			errpatch[npatch][0] = errormin[0];
			int nn = 1;
			for (int i = 0; i < MAX - 1 && nn < CONTENTS; i++)
			{
				bool ok = true;
				for (int j = 0; j < nn; j++)
					if (idscale[npatch][j] == smin[i] && idrot[npatch][j] == rotmin[i]
						&& (int)sqrt((double)((deltax[npatch][j] - dminx[i])*(deltax[npatch][j] - dminx[i]) + (deltay[npatch][j] - dminy[i])*(deltay[npatch][j] - dminy[i]))) < RADIUS) ok = false;
				if (ok)
				{
					deltax[npatch][nn] = dminx[i]; deltay[npatch][nn] = dminy[i];
					scalex[npatch][nn] = sfactorx[smin[i]]; scaley[npatch][nn] = sfactory[smin[i]];
					idscale[npatch][nn] = smin[i];
					idrot[npatch][nn] = rotmin[i];
					errpatch[npatch][nn] = errormin[i];
					nn++;
				}
			}

			if (nn < CONTENTS) {
				for (int i = nn; i < CONTENTS; i++)
				{ deltax[npatch][i] = deltax[npatch][i%nn]; deltay[npatch][i] = deltay[npatch][i%nn]; }
			}

			for (int i = 0; i < CONTENTS; i++) 
			{
				dminx[i] = deltax[npatch][i]; dminy[i] = deltay[npatch][i];
				smin[i] = idscale[npatch][i];
				rotmin[i] = idrot[npatch][i];
			}
		}

		// refining translation vectors to resolution of tile
		////////////////////////////////////////////////////
		for (int npatch = 0; npatch < NPATCH; npatch++)
		{
			hvBitmap patchb(2 * tile.sizeX(), 2 * tile.sizeY(), false);
			hvBitmap feature(2 * tile.sizeX(), 2 * tile.sizeY(), false);
			for (int j = 0; j < patchfrags[npatch]->size(); j++)
			{
				int k;
				int fragind = patchfrags[npatch]->at(j);
				for (k = 0; k < f.fragmentSize(fragind); k++)
				{
					int u, v; f.pixelFragment(fragind, k, u, v);
					hvVec4<int> col = patchpf.get(u, v);
					int x = u + (col.W() == 3 || col.W() == 1 ? tile.sizeX() : 0);
					int y = v + (col.W() == 3 || col.W() == 2 ? tile.sizeY() : 0);
					patchb.set(x, y, true);
				}
			}
			patchb.fillholes(); feature = patchb;
			feature.erosion(3, 3); feature.erosion(3, 3); ~feature;
			patchb.dilatation(3, 3); patchb.dilatation(3, 3); patchb.dilatation(3, 3);
			//~feature; feature.dilatation(3,3); feature.fillholes(); ~feature;
			patchb &= feature;
			std::vector<hvVec2<int> > border(100000); border.clear();
			int minx = 2 * tile.sizeX(), miny = 2 * tile.sizeY(), maxx = 0, maxy = 0;
			for (int i = 0; i < 2 * tile.sizeX(); i++)
			{
				for (int j = 0; j < 2 * tile.sizeY(); j++)
				{
					if (patchb.get(i, j))
					{
						border.push_back(hvVec2<int>(i, j));
						if (minx > i) minx = i;
						if (maxx < i) maxx = i;
						if (miny > j) miny = j;
						if (maxy < j) maxy = j;
					}
				}
			}
			ccx[npatch] = (maxx + minx) / 2; ccy[npatch] = (maxy + miny) / 2;

			if (border.size() > 100000) std::cerr << "too many border pixels!\n" << std::endl;
			//int *bordneigh = new int[100000 * 50];
			//int *nbordneigh = new int[100000];
			std::vector<int> bordneigh(100000 * 50);
			std::vector<int> nbordneigh(100000);

			for (int i = 0; i < border.size(); i++) nbordneigh[i] = 0;
			for (int i = 0; i < border.size(); i++) for (int j = 0; j < border.size(); j++)
			{
				hvVec2<int> posa = border.at(i);
				hvVec2<int> posb = border.at(j);
				if (sqrt((double)((posa.X() - posb.X())*(posa.X() - posb.X()) + (posa.Y() - posb.Y())*(posa.Y() - posb.Y()))) < 2)
				{
					if (nbordneigh[i] < 50)
					{
						bordneigh[i * 50 + nbordneigh[i]] = j;
						nbordneigh[i]++;
					}
					else { hvFatal("too many neighbors!\n"); }
				}
			}
			
			for (int nexg = 0; nexg < CONTENTS; nexg++)
			{
				int dminx, dminy;
				double errormin = 1000000000.0;
				for (int i = 0; i < ratio + 4; i += 2) for (int j = 0; j < ratio + 4; j += 2)
				{
					int ii = deltax[npatch][nexg] * ratio + i - 2;
					int jj = deltay[npatch][nexg] * ratio + j - 2;
					int count = 0;
					double errsum = 0.0;
					for (int k = 0; k < border.size(); k++)
					{
						hvVec2<int> pos = border.at(k);
						double theta = rotangles[idrot[npatch][nexg]];
						int bx = (int)((double)(pos.X() - ccx[npatch])*scalex[npatch][nexg] * cos(theta) - (double)(pos.Y() - ccy[npatch])*scaley[npatch][nexg] * sin(theta)) + ii;
						int by = (int)((double)(pos.X() - ccx[npatch])*scalex[npatch][nexg] * sin(theta) + (double)(pos.Y() - ccy[npatch])*scaley[npatch][nexg] * cos(theta)) + jj;
						while (bx < 0) bx += tile.sizeX(); while (bx >= tile.sizeX()) bx -= tile.sizeX();
						while (by < 0) by += tile.sizeY(); while (by >= tile.sizeY()) by -= tile.sizeY();
						int ppx = pos.X(); ppx = ppx%tile.sizeX();
						int ppy = pos.Y(); ppy = ppy%tile.sizeY();
						hvColRGB<unsigned char> cola = tile.get(bx, by);
						hvColRGB<unsigned char> colb = tile.get(ppx, ppy);
						//if (bx<1||bx>=example.sizeX()-1||by<1||by>=example.sizeY()-1) ok=false;
						//else
						//errsum+=cola.squaredDifferenceNorm(255.0,colb);

						if (!(tilesal.get(ppx, ppy).equals(tilesal.get(bx, by))))
						{
							count++;
							double err = 0.0;
							for (int q = 0; q < nbordneigh[k]; q++)
							{
								hvVec2<int> qpos = border.at(bordneigh[k * 50 + q]);
								double theta = rotangles[idrot[npatch][nexg]];
								int qbx = (int)((double)(qpos.X() - ccx[npatch])*scalex[npatch][nexg] * cos(theta) - (double)(qpos.Y() - ccy[npatch])*scaley[npatch][nexg] * sin(theta)) + ii;
								int qby = (int)((double)(qpos.X() - ccx[npatch])*scalex[npatch][nexg] * sin(theta) + (double)(qpos.Y() - ccy[npatch])*scaley[npatch][nexg] * cos(theta)) + jj;
								while (qbx < 0) qbx += tile.sizeX(); while (qbx >= tile.sizeX()) qbx -= tile.sizeX();
								while (qby < 0) qby += tile.sizeY(); while (qby >= tile.sizeY()) qby -= tile.sizeY();
								//if (qbx<1||qbx>=example.sizeX()-1||qby<1||qby>=example.sizeY()-1) ok=false;
								//else
								{
									int ppx = qpos.X(); ppx = ppx%tile.sizeX();
									int ppy = qpos.Y(); ppy = ppy%tile.sizeY();
									double mm = tilep.get(ppx, ppy).squaredDifferenceNorm(255.0, tilep.get(qbx, qby));
									err += mm;
								}
							}
							//errsum+=err/(double)nbordneigh[k];
							errsum += err;
						}

					}
					if (errormin > errsum)
					{
						errormin = errsum; dminx = ii; dminy = jj;
					}
				}
				deltax[npatch][nexg] = dminx;
				deltay[npatch][nexg] = dminy;
				errpatch[npatch][nexg] = errormin;
			}
		}

		for (int npatch = 0; npatch < NPATCH; npatch++)
		{
			for (int i = 0; i < CONTENTS; i++)
			{
				for (int j = i + 1; j < CONTENTS; j++)
				{
					if (errpatch[npatch][j] < errpatch[npatch][i])
					{
						int dminx = deltax[npatch][i], dminy = deltay[npatch][i], smin = idscale[npatch][i], rotmin = idrot[npatch][i];
						double errormin = errpatch[npatch][i];
						deltax[npatch][i] = deltax[npatch][j]; deltay[npatch][i] = deltay[npatch][j];
						idscale[npatch][i] = idscale[npatch][j];
						idrot[npatch][i] = idrot[npatch][j];
						errpatch[npatch][i] = errpatch[npatch][j];
						deltax[npatch][j] = dminx; deltay[npatch][j] = dminy;  idscale[npatch][j] = smin; idrot[npatch][j] = rotmin;
						errpatch[npatch][j] = errormin;
					}
				}
			}

		}

		///////////////////////////////////////////////////////////////////////////////////////////

		auto compute_turbul = [&](double turb_scale, double temps, int PERIOD, int FSCALE, const std::string& ext)
		{
			hvPictRGB<unsigned char> nouv(tile.sizeX()*PERIOD / FSCALE, tile.sizeY()*PERIOD / FSCALE, hvColRGB<unsigned char>());
			hvPictRGB<unsigned char> nouvorig(tile.sizeX()*PERIOD / FSCALE, tile.sizeY()*PERIOD / FSCALE, hvColRGB<unsigned char>());
			hvPictRGB<unsigned char> nouvborder(tile.sizeX()*PERIOD / FSCALE, tile.sizeY()*PERIOD / FSCALE, hvColRGB<unsigned char>());

			int nbj = tile.sizeY()*PERIOD / FSCALE;
			int nbi = tile.sizeX()*PERIOD / FSCALE;
			for (int ii = 0; ii < nbi; ii++)
			{
				for (int jj = 0; jj < nbj; jj++)
				{
					int oi = tile.sizeX() + ii*FSCALE;
					int oj = tile.sizeY() + jj*FSCALE;

					double ttx = proc.turbulence(turb_scale*(double)oi, turb_scale*(double)oj);
					double tty = proc.turbulence(turb_scale*(double)oi + 4.543, turb_scale*(double)oj + 9.654);
					int i = tile.sizeX() + ii*FSCALE + (int)(ttx*TAMPX*temps);
					int j = tile.sizeY() + jj*FSCALE + (int)(tty*TAMPY*temps);

					hvColRGB<unsigned char> col;
					hvColRGB<double> cc;
					std::vector<hvColRGB<double>> coln(NCLUSTERS);
					hvVec4<int> vv = patchpf.get(i%tile.sizeX(), j%tile.sizeY());
					int npatch = vv.X();
					int pcenterx = (i / tile.sizeX())*tile.sizeX() + ccx[npatch] - (vv.W() == 3 || vv.W() == 1 ? tile.sizeX() : 0);
					int pcentery = (j / tile.sizeY())*tile.sizeY() + ccy[npatch] - (vv.W() == 3 || vv.W() == 2 ? tile.sizeY() : 0);

					hvNoise::seeding(pcenterx / ratio / 2, pcentery / ratio / 2, 10);
					int ind = (int)(0.5f*(hvNoise::next() + 1.0f)*(float)(CONTENTS));
					ind = ind%CONTENTS;

					int pi, pj;
					double theta = rotangles[idrot[npatch][ind]];
					pi = (int)((double)(i - pcenterx)*sfactorx[idscale[npatch][ind]] * cos(theta) - (double)(j - pcentery)*sfactory[idscale[npatch][ind]] * sin(theta)) + deltax[npatch][ind]; while (pi < 0) pi = pi + tile.sizeX(); while (pi >= tile.sizeX()) pi -= tile.sizeX();
					pj = (int)((double)(i - pcenterx)*sfactorx[idscale[npatch][ind]] * sin(theta) + (double)(j - pcentery)*sfactory[idscale[npatch][ind]] * cos(theta)) + deltay[npatch][ind]; while (pj < 0) pj = pj + tile.sizeY(); while (pj >= tile.sizeY()) pj -= tile.sizeY();

					col = tile.get(pi, pj);
					hvColRGB<unsigned char> coltp = tile.get(i%tile.sizeX(), j%tile.sizeY());
					double alpha = (double)patcherod.get(i%tile.sizeX(), j%tile.sizeY()).RED() / 255.0;
					double ww = 1.0 - (double)weight.get(pi, pj).RED() / 255.0;
					hvColRGB<unsigned char> colint;
					if (sqrt(coltp.squaredDifference(col) / 3.0) < COL_THRESH) colint.blend(col, coltp, pow(alpha, ww)); else colint = col;
					col = colint;
					cc = (hvColRGB<double>)col;
					for (int ind = 0; ind < NCLUSTERS; ind++)
					{
						int nx = oi%noisei[ind].sizeX();
						int ny = oj%noisei[ind].sizeY();
						coln[ind] = (hvColRGB<double>)noisei[ind].get(nx, ny);

						coln[ind] -= hvColRGB<double>(127.0);
						coln[ind].scale((double)noisetile[ind].get(pi, pj).RED() / 255.0);
						cc += coln[ind];
					}
					cc.clamp(0.0, 255.0);
					nouv.update(ii, jj, (hvColRGB<unsigned char>)cc);

					/////////////////////////////////////
					col = tilep.get(pi, pj);
					coltp = tilep.get(i%tile.sizeX(), j%tile.sizeY());
					if (sqrt(coltp.squaredDifference(col) / 3.0) < COL_THRESH) colint.blend(col, coltp, pow(alpha, ww)); else colint = col;
					col = colint;
					cc = (hvColRGB<double>)col;
					double wwp = 0.0;
					hvColRGB<double> ccn(0.0);
					for (int ind = 0; ind < NCLUSTERS; ind++)
					{
						int nx = oi%noisei[ind].sizeX();
						int ny = oj%noisei[ind].sizeY();
						coln[ind] = (hvColRGB<double>)noisei[ind].get(nx, ny);
						coln[ind] -= hvColRGB<double>(127.0);
						wwp += (double)noisetileorig[ind].get(pi, pj) / 255.0;
						coln[ind].scale((double)noisetileorig[ind].get(pi, pj) / 255.0);
						ccn += coln[ind];
					}
					ccn.scale(1.0 / wwp);
					cc += ccn;
					cc.clamp(0.0, 255.0);
					nouvborder.update(ii, jj, (hvColRGB<unsigned char>)cc);
					hvColRGB<unsigned char> colt = tilep.get(pi, pj);
					cc = (hvColRGB<double>)colt;
					cc.clamp(0.0, 255.0);
					colt = (hvColRGB<unsigned char>)cc;
					cc.difference(nouv.get(ii, jj), colt, 127.0, 127.0);
					colt = tileorig.get(pi, pj);
					cc = (hvColRGB<double>)colt;
					cc.clamp(0.0, 255.0);
					colt = (hvColRGB<unsigned char>)cc;
					nouvorig.update(ii, jj, colt);
				}
			}

			std::string name = base_dir + name_file + "_Final_blt_boundary_" + ext + ".png";
			ASTex::save_from_hv(name, nouv);

			name = base_dir + name_file + "_Final_blt_" + ext + ".png";
			ASTex::save_from_hv(name, nouvborder);

			name = base_dir + name_file + "_Final_no_blt" + ext + ".png";
			ASTex::save_from_hv(name, nouvorig);
		};


		std::thread thr1([&]() { compute_turbul(0.0, 0.0, SZ_MULT, 1, ""); });
		std::thread thr2([&]() { compute_turbul(0.5, 0.5, SZ_MULT, 1, "_turb"); });

		thr1.join();
		thr2.join();

}

#if defined(__clang__)
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#elif defined(_MSC_VER)
#pragma warning( pop )
#endif
