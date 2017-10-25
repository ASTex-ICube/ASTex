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



#include <cmath>
#include <ASTex/special_io.h>
#include <ASTex/fourier.h>
#include <ASTex/local_spectrum.h>
#include <ASTex/utils.h>
#include <ASTex/colorspace_filters.h>
#include <ASTex/mask.h>
#include <ASTex/distances_maps.h>
#include <ASTex/pca.h>
#include <ctime>


namespace ASTex
{

template<typename MASK>
int true_prop_verif(const ImageRGBd& /*input*/, int acc_size, int block_x, int block_y, const MASK& mask)
{
	const int T = acc_size;

	// compute auto-correlation
	long int mask_size = 0;
	for (int y = 0; y < T; ++y)
	{
		for (int x = 0; x < T; ++x)
		{
			if (mask(block_x + x, block_y + y))
			{
				++mask_size;
			}
		}
	}
	return mask_size;
}

void computeTransfer(std::vector<uint8_t>& trans, std::vector<int>& histo, std::vector<int>& proc_histo, int i, int j, int sumh, int a, int b, int sumr)
{

	if (i == j)
	{
		for (int k = a; k <= b; k++)
			trans[k] = (uint8_t)i;
		return;
	}

	int ind=0;
	if (a == b)
	{
		for (ind = i; histo[ind] == 0; ind++);

		for (int k = a; k <= b; k++)
			trans[k] = (uint8_t)ind;

		return;
	}

	int count = 0;

	for (int k = i; k <= j; k++)
		if (histo[k] != 0)
		{
			count++;
			ind = k;
		}

	if (count == 1)
	{
		for (int k = a; k <= b; k++)
			trans[k] = (uint8_t)ind;
		return;
	}

	count = 0;
	for (int k = a; k <= b; k++)
		if (proc_histo[k] != 0)
		{
			count++;
			ind = k;
		}

	if (count == 1)
	{
		for (int k = a; k <= b; k++)
			trans[k] = (uint8_t)((i + j) / 2);
		return;
	}

	int k;
	int sumnh = 0;
	int sumnr = 0;

	for (k = i; k <= j && sumnh < sumh / 2; k++)
	{
		sumnh += histo[k];
	}

	int ii = k - 1;
	if (sumh - sumnh == 0)
	{
		sumnh -= histo[ii];
		ii--;
	}

	for (k = a; k <= b && sumnr < sumr / 2; k++)
	{
		sumnr += proc_histo[k];
	}

	int aa = k - 1; if (sumr - sumnr == 0) { sumnr -= proc_histo[aa]; aa--; }


	computeTransfer(trans, histo, proc_histo, i, ii, sumnh, a, aa, sumnr);
	computeTransfer(trans, histo, proc_histo, ii + 1, j, sumh - sumnh, aa + 1, b, sumr - sumnr);
}


void Synthesis(const std::string& out_path, const std::string& inputfile,/* const std::string& synthese_name,*/
	/*const ImageRGBd& N_input,*/ const ImageRGBd& S_input,
	const std::vector<PCA>& pcas,
	const std::vector<std::vector<ImageSpectrald>>& spectrum,
	const std::vector<ImageGrayd>& guidance, std::vector<std::vector<ImageGrayd>> coords_per_mask)
{


	std::string filename_plus_ext = ASTex::IO::remove_path(inputfile);
	std::string name_file = ASTex::IO::remove_ext(filename_plus_ext);

	ImageRGBd input;
	input.load(inputfile);

	typedef ImageRGBd::ItkImg IMG_DBL;
	// first transform [0,255] double -> [0,1] double
	ColorSpace::FilterRGB255To01<IMG_DBL, IMG_DBL>::Pointer filter0 =
		ColorSpace::FilterRGB255To01<IMG_DBL, IMG_DBL>::New();
	filter0->SetInput(input.itk());
	filter0->Update();
	input.itk() = filter0->GetOutput();



	ImageRGBd noise_localised;
	noise_localised.initItk(S_input.width(), S_input.height(), true);

	ImageRGBd input_synthesis;
	input_synthesis.initItk(S_input.width(), S_input.height(), true);

	int GEOFFREY_MAGIC_SIZE = 8 * std::max(S_input.width(), S_input.height()); //4000 ?

	std::vector<std::vector<ImageGrayd>> noise_synthesis_per_cord_per_mask;
	noise_synthesis_per_cord_per_mask.resize(3, std::vector<ImageGrayd>(guidance.size()));
	for (uint32_t i = 0; i < 3; ++i)
		for (uint32_t j = 0; j < guidance.size(); ++j)
			noise_synthesis_per_cord_per_mask[i][j].initItk(GEOFFREY_MAGIC_SIZE, GEOFFREY_MAGIC_SIZE);
	//    noise_synthesis_per_cord_per_mask[i][j].initItk(S_input.width(),S_input.height());


	std::vector<ImageRGBd> noise_synthesis;
	noise_synthesis.resize(guidance.size());


	for (uint32_t i = 0; i < guidance.size(); ++i)
	{
		//        noise_synthesis[i].initItk(S_input.width(),S_input.height(),true);
		noise_synthesis[i].initItk(GEOFFREY_MAGIC_SIZE, GEOFFREY_MAGIC_SIZE, true);
	}

	srand(time(NULL));

	for (uint32_t i = 0; i < guidance.size(); ++i)
	{
		Fourier::RPnoise_mosaic_bandes(spectrum[0][i], noise_synthesis_per_cord_per_mask[0][i]);
		Fourier::RPnoise_mosaic_bandes(spectrum[1][i], noise_synthesis_per_cord_per_mask[1][i]);
		Fourier::RPnoise_mosaic_bandes(spectrum[2][i], noise_synthesis_per_cord_per_mask[2][i]);

		// EGALISATION D'HISTOGRAMME

		//        int lol = coords_per_mask[0][i].width()*coords_per_mask[0][i].height();

		int taille_h = noise_synthesis_per_cord_per_mask[0][i].height();
		int taille_w = noise_synthesis_per_cord_per_mask[0][i].width();

		//		int histovar[256], shistovar[512];
		std::vector<int32_t> histovar(256, 0);
		std::vector<int32_t> shistovar(512, 0);
		std::vector<uint8_t> transvar(512, 0);

		//		for (int x = 0; x < 256; x++)
		//			histovar[x] = 0;

		//		for (int x = 0; x < 512; x++)
		//		{
		//			shistovar[x] = 0;
		//			transvar[x]= 0;
		//		}

		float Imin = 1; float Imax = 0;
		float smin = 1; float smax = 0;
		//min max input

		for (int x = 0; x < coords_per_mask[0][i].width(); x++) for (int y = 0; y < coords_per_mask[0][i].height(); y++)
		{
			if (guidance[i].pixelAbsolute(x, y) == 1)
			{
				if (coords_per_mask[0][i].pixelAbsolute(x, y) < Imin) Imin = coords_per_mask[0][i].pixelAbsolute(x, y);
				if (coords_per_mask[0][i].pixelAbsolute(x, y) > Imax) Imax = coords_per_mask[0][i].pixelAbsolute(x, y);
			}

			if (noise_synthesis_per_cord_per_mask[0][i].pixelAbsolute(x, y) < smin) smin = noise_synthesis_per_cord_per_mask[0][i].pixelAbsolute(x, y);
			if (noise_synthesis_per_cord_per_mask[0][i].pixelAbsolute(x, y) > smax) smax = noise_synthesis_per_cord_per_mask[0][i].pixelAbsolute(x, y);

		}


		for (int x = 0; x < coords_per_mask[0][i].width(); x++) for (int y = 0; y < coords_per_mask[0][i].height(); y++)
		{
			if (guidance[i].pixelAbsolute(x, y) == 1)
			{
				//COMPOSANTE 0
				int bin = (int)((coords_per_mask[0][i].pixelAbsolute(x, y) - Imin) / (Imax - Imin)*256.0);
				if (bin < 0) bin = 0; else if (bin >= 256) bin = 255; 
				histovar[bin] += 1;
			}
		}

		for (int x = 0; x < noise_synthesis_per_cord_per_mask[0][i].width(); x++) for (int y = 0; y < noise_synthesis_per_cord_per_mask[0][i].height(); y++)
		{

			int bin = (int)((noise_synthesis_per_cord_per_mask[0][i].pixelAbsolute(x, y) - smin) / (smax - smin)*512.0);
			if (bin < 0) bin = 0; else if (bin >= 512) bin = 511;
			shistovar[bin] += 1;


		}
		int hcount = 0;

		for (int k = 0; k < 256; k++)
			hcount += histovar[k];

		computeTransfer(transvar, histovar, shistovar, 0, 255, hcount, 0, 511, taille_w*taille_h);

		for (int y = 0; y < taille_h; y++)
			for (int x = 0; x < taille_w; x++)
			{

				double vvar = (noise_synthesis_per_cord_per_mask[0][i].pixelAbsolute(x, y) - smin) / (smax - smin)*512.0;

				if (vvar < 0.0)
					vvar = 0.0;
				else
					if (vvar >= 511.0)
						vvar = 511.0;

				vvar = ((double)transvar[(int)vvar] / 255.0*(double)(Imax - Imin) + (double)Imin);

				noise_synthesis_per_cord_per_mask[0][i].pixelAbsolute(x, y) = vvar;

			}

		///////////////////////////////////////////////////////GREEN


//		for (int x = 0; x < 256; x++)
//			histovar[x] = 0;

//		for (int x = 0; x < 512; x++)
//		{
//			shistovar[x] = 0;
//			transvar[x] = 0;
//		}
		histovar.assign(256, 0);
		shistovar.assign(512, 0);
		transvar.assign(512, 0);


		Imin = 1;  Imax = 0;
		smin = 1;  smax = 0;
		//min max input
		for (int x = 0; x < coords_per_mask[1][i].width(); x++) for (int y = 0; y < coords_per_mask[1][i].height(); y++)
		{
			if (guidance[i].pixelAbsolute(x, y) == 1)
			{
				if (coords_per_mask[1][i].pixelAbsolute(x, y) < Imin) Imin = coords_per_mask[1][i].pixelAbsolute(x, y);
				if (coords_per_mask[1][i].pixelAbsolute(x, y) > Imax) Imax = coords_per_mask[1][i].pixelAbsolute(x, y);
			}


			if (noise_synthesis_per_cord_per_mask[1][i].pixelAbsolute(x, y) < smin) smin = noise_synthesis_per_cord_per_mask[1][i].pixelAbsolute(x, y);
			if (noise_synthesis_per_cord_per_mask[1][i].pixelAbsolute(x, y) > smax) smax = noise_synthesis_per_cord_per_mask[1][i].pixelAbsolute(x, y);
		}


		for (int x = 0; x < coords_per_mask[1][i].width(); x++) for (int y = 0; y < coords_per_mask[1][i].height(); y++)
		{
			if (guidance[i].pixelAbsolute(x, y) == 1)
			{
				//COMPOSANTE 0
				int bin = (int)((coords_per_mask[1][i].pixelAbsolute(x, y) - Imin) / (Imax - Imin)*256.0);
				if (bin < 0) bin = 0; else if (bin >= 256) bin = 255; 
				histovar[bin] += 1;
			}
		}

		for (int x = 0; x < noise_synthesis_per_cord_per_mask[1][i].width(); x++) for (int y = 0; y < noise_synthesis_per_cord_per_mask[1][i].height(); y++)
		{

			int bin = (int)((noise_synthesis_per_cord_per_mask[1][i].pixelAbsolute(x, y) - smin) / (smax - smin)*512.0);
			if (bin < 0) bin = 0; else if (bin >= 512) bin = 511;
			shistovar[bin] += 1;


		}
		hcount = 0;

		for (int k = 0; k < 256; k++) hcount += histovar[k];

		computeTransfer(transvar, histovar, shistovar, 0, 255, hcount, 0, 511, taille_w*taille_h);

		for (int x = 0; x < taille_w; x++) for (int y = 0; y < taille_h; y++)
		{

			double vvar = (noise_synthesis_per_cord_per_mask[1][i].pixelAbsolute(x, y) - smin) / (smax - smin)*512.0;

			if (vvar < 0.0)
				vvar = 0.0;
			else if (vvar >= 511.0)
				vvar = 511.0;

			vvar = ((double)transvar[(int)vvar] / 255.0*(double)(Imax - Imin) + (double)Imin);

			noise_synthesis_per_cord_per_mask[1][i].pixelAbsolute(x, y) = vvar;

		}

		/////////////////////////////////////////////////////////BLUE
		///
		///
//		for (int x = 0; x < 256; x++)
//			histovar[x] = 0;
//		for (int x = 0; x < 512; x++)
//		{
//			shistovar[x] = 0;
//			transvar[x]=0;
//		}

		histovar.assign(256, 0);
		shistovar.assign(512, 0);
		transvar.assign(512, 0);


		Imin = 1;  Imax = 0;
		smin = 1;  smax = 0;
		//min max input
		for (int y = 0; y < coords_per_mask[2][i].height(); y++)
			for (int x = 0; x < coords_per_mask[2][i].width(); x++)
			{
				if (guidance[i].pixelAbsolute(x, y) == 1)
				{
					if (coords_per_mask[2][i].pixelAbsolute(x, y) < Imin) Imin = coords_per_mask[2][i].pixelAbsolute(x, y);
					if (coords_per_mask[2][i].pixelAbsolute(x, y) > Imax) Imax = coords_per_mask[2][i].pixelAbsolute(x, y);

				}

				if (noise_synthesis_per_cord_per_mask[2][i].pixelAbsolute(x, y) < smin) smin = noise_synthesis_per_cord_per_mask[2][i].pixelAbsolute(x, y);
				if (noise_synthesis_per_cord_per_mask[2][i].pixelAbsolute(x, y) > smax) smax = noise_synthesis_per_cord_per_mask[2][i].pixelAbsolute(x, y);
			}


		for (int y = 0; y < coords_per_mask[2][i].height(); y++)
			for (int x = 0; x < coords_per_mask[2][i].width(); x++)
			{
				if (guidance[i].pixelAbsolute(x, y) == 1)
				{
					//COMPOSANTE 0
					int bin = (int)((coords_per_mask[2][i].pixelAbsolute(x, y) - Imin) / (Imax - Imin)*256.0);
					if (bin < 0) bin = 0; else if (bin >= 256) bin = 255;
					histovar[bin] += 1;
				}
			}

		for (int y = 0; y < noise_synthesis_per_cord_per_mask[2][i].height(); y++)
			for (int x = 0; x < noise_synthesis_per_cord_per_mask[2][i].width(); x++)
			{

				int bin = (int)((noise_synthesis_per_cord_per_mask[2][i].pixelAbsolute(x, y) - smin) / (smax - smin)*512.0);
				if (bin < 0) bin = 0; else if (bin >= 512) bin = 511;
				shistovar[bin] += 1;


			}
		hcount = 0;

		for (int k = 0; k < 256; k++)
			hcount += histovar[k];

		computeTransfer(transvar, histovar, shistovar, 0, 255, hcount, 0, 511, taille_w*taille_h);

		for (int y = 0; y < taille_h; y++)
			for (int x = 0; x < taille_w; x++)
			{

				double vvar = (noise_synthesis_per_cord_per_mask[2][i].pixelAbsolute(x, y) - smin) / (smax - smin)*512.0;

				if (vvar < 0.0)
					vvar = 0.0;
				else if (vvar >= 511.0)
					vvar = 511.0;

				vvar = ((double)transvar[(int)vvar] / 255.0*(double)(Imax - Imin) + (double)Imin);

				noise_synthesis_per_cord_per_mask[2][i].pixelAbsolute(x, y) = vvar;

			}

		// PROJECTION INVERSE DANS NOISE_SYNTHESIS
		pcas[i].back_project(noise_synthesis_per_cord_per_mask[0][i], noise_synthesis_per_cord_per_mask[1][i], noise_synthesis_per_cord_per_mask[2][i], noise_synthesis[i]);


		IO::save(noise_synthesis[i], out_path + name_file + "_noise_synth" + std::to_string(i) + ".png", 1.0, 0.5, 0.5, 0.5);
	}

}



void Synthesis_corel(const std::string& filename_source, const std::string& base_dir, int nb_clusters)
{
	std::string name_file = IO::remove_ext(IO::remove_path(filename_source));
	std::string in_noise = base_dir + name_file + "_filtered_N.png";
	std::string in_struct = base_dir + name_file + "_filtered_S.png";
	std::vector<std::string> dist_maps(nb_clusters);
	std::vector<std::string> guidance_maps(nb_clusters);
	for (int i = 0; i < nb_clusters; ++i)
	{
		dist_maps[i] = base_dir + name_file + "_mask_dist" + std::to_string(i) + ".png";
		guidance_maps[i] = base_dir + name_file + "_mask_bin" + std::to_string(i) + ".png";
	}

	assert(dist_maps.size() == guidance_maps.size());
	int nb_maps = dist_maps.size();

	ImageRGBd input;
	input.load(filename_source);

	ImageRGBd N_input;
	N_input.load(in_noise);

	ImageRGBd S_input;
	S_input.load(in_struct);


	typedef ImageRGBd::ItkImg IMG_DBL;
	// first transform [0,255] double -> [0,1] double
	ColorSpace::FilterRGB255To01<IMG_DBL, IMG_DBL>::Pointer filter0 =
		ColorSpace::FilterRGB255To01<IMG_DBL, IMG_DBL>::New();

	ColorSpace::FilterRGB255To01<IMG_DBL, IMG_DBL>::Pointer filter1 =
		ColorSpace::FilterRGB255To01<IMG_DBL, IMG_DBL>::New();

	ColorSpace::FilterRGB255To01<IMG_DBL, IMG_DBL>::Pointer filter2 =
		ColorSpace::FilterRGB255To01<IMG_DBL, IMG_DBL>::New();

	filter0->SetInput(input.itk());
	filter0->Update();
	input.itk() = filter0->GetOutput();

	filter1->SetInput(N_input.itk());
	filter1->Update();
	N_input.itk() = filter1->GetOutput();

	filter2->SetInput(S_input.itk());
	filter2->Update();
	S_input.itk() = filter2->GetOutput();

	for (int i = 0; i < N_input.width(); ++i)
		for (int j = 0; j < N_input.height(); ++j)
		{

			N_input.pixelAbsolute(i, j)[0] -= 0.5;
			N_input.pixelAbsolute(i, j)[1] -= 0.5;
			N_input.pixelAbsolute(i, j)[2] -= 0.5;
		}

	std::vector<std::vector<ImageSpectrald>> spectrum;
	spectrum.resize(3, std::vector<ImageSpectrald>(nb_maps));

	std::vector<ImageGrayd> guidance;
	guidance.resize(nb_maps);

	std::vector<ImageGrayd> guidance_synthesis;
	guidance_synthesis.resize(nb_maps);


	for (int i = 0; i < nb_maps; ++i)
	{
		IO::load(guidance[i], dist_maps[i]);
		IO::load(guidance_synthesis[i], guidance_maps[i]);
	}


	std::vector<int> nb_pix_mask;
	nb_pix_mask.resize(nb_maps);

	std::vector<double> prop_selection_masque;
	prop_selection_masque.resize(nb_maps);


	for (uint32_t k = 0; k < guidance.size(); ++k) {
		prop_selection_masque[k] = 0;
		nb_pix_mask[k] = 0;
	}

	for (uint32_t k = 0; k < guidance.size(); ++k) {
		for (int i = 0; i < N_input.width(); ++i)
			for (int j = 0; j < N_input.height(); ++j)
			{
				if (guidance_synthesis[k].pixelAbsolute(i, j) == 1)
					nb_pix_mask[k]++;
			}
	}
	int nb_pix_max = N_input.width()*N_input.height();
	for (uint32_t k = 0; k < guidance.size(); ++k) {
		prop_selection_masque[k] = double(nb_pix_mask[k]) / double(nb_pix_max);
	}

	//ACP

	//OUTPUT PCA

	std::vector<std::vector<ImageGrayd>> coords_per_mask;
	coords_per_mask.resize(3, std::vector<ImageGrayd>(guidance.size()));
	for (int i = 0; i < 3; ++i)
		for (uint32_t j = 0; j < guidance.size(); ++j)
			coords_per_mask[i][j].initItk(input.width(), input.height(), true);

	/********************* PCA ************************/

	std::vector<PCA> pcas;
	for (uint32_t k = 0; k < guidance.size(); ++k)
		pcas.push_back(PCA(N_input)); // Input : sensé être notre noise RGB

	std::vector<ImageGrayd> masks_normalized;
	masks_normalized.resize(guidance.size());

	normalize_distance_maps(guidance, masks_normalized);

	std::vector<MaskSmallestValue<double>> masks_binar;
	//    masks_binar.resize(guidance.size());


	ImageGrayd inputGrey;
	IO::load_RGB_2_gray(inputGrey, filename_source);

	//Pour chaque masque
	//    int AC_size =64;
	for (uint32_t k = 0; k < guidance.size(); ++k)
	{
		//Creation de notre masque et export du binaire
		float val = prop_selection_masque[k] - (0.1*prop_selection_masque[k]);

		MaskSmallestValue<double> m(guidance[k], val);
		masks_binar.push_back(m);


		//Calcul de la base
		pcas[k].computePCA(m);

		pcas[k].project(coords_per_mask[0][k], coords_per_mask[1][k], coords_per_mask[2][k]);

		ImageGrayu8 mask_tmp;
		mask_tmp.initItk(input.width(), input.height(), true);
		m.export_binary_image(mask_tmp);
		mask_tmp.save(base_dir + "BINARY_nbclust_" + std::to_string(k) + ".png");

	}


	std::vector<int>  AC_sizes;
	std::vector<double>  seuil;
	std::vector<bool>  seuil_search_canal;

	for (std::size_t k = 0; k < guidance.size(); ++k)
	{
		seuil.push_back(0.9);
		seuil_search_canal.push_back(false);
		AC_sizes.push_back(64);
	}


	// TODO SEUIL POUR CHAQUE ZONE PASSEZ EN VEC


	bool seuil_search = true;
	int sum_crop[3];

	while (seuil_search) {
		seuil_search = false;

		for (uint32_t k = 0; k < guidance.size(); ++k) {
			if (seuil_search_canal[k] == false) {
				sum_crop[k] = 0;
				int nb_crop = 0;
				for (int block_y = 1; block_y + AC_sizes[k] <= input.height(); block_y += 1)
				{
					for (int block_x = 1; block_x + AC_sizes[k] <= input.width(); block_x += 1)
					{
						nb_crop++;
						int mask_size = true_prop_verif(input, AC_sizes[k], block_x, block_y, masks_binar[k]);

						float prop = mask_size / double(AC_sizes[k] * AC_sizes[k]);
						if (prop > seuil[k])
							sum_crop[k]++;
					}
				}
				//        sum_crop/=nb_crop;
			}
		}
		for (uint32_t k = 0; k < guidance.size(); ++k) {

			if (sum_crop[k] < AC_sizes[k] * AC_sizes[k] || sum_crop[k] < 1500)
			{
				seuil[k] -= 0.1;

				if (seuil[k] < 0.2 && AC_sizes[k] == 8) {
					seuil_search_canal[k] = true;
				}
				else if (seuil[k] < 0.1) {
					seuil[k] = 0.9;
					AC_sizes[k] /= 2;
				}
				else {
					seuil_search = true;
				}
			}
			else
			{
				seuil_search_canal[k] = true;
			}
		}
		for (uint32_t k = 0; k < guidance.size(); ++k) {
			if (seuil_search_canal[k] == false)
				seuil_search = true;
		}

	}

	for (uint32_t k = 0; k < guidance.size(); ++k)
	{
		// Charger les spectres par auto-covariance
		spectrum[1][k].initItk(AC_sizes[k], AC_sizes[k], true);
		spectrum[0][k].initItk(AC_sizes[k], AC_sizes[k], true);
		spectrum[2][k].initItk(AC_sizes[k], AC_sizes[k], true);
	}


	for (uint32_t k = 0; k < guidance.size(); ++k)
	{
		Fourier::spectrum_by_autocorrelation_small_size(coords_per_mask[0][k], spectrum[0][k], masks_binar[k], seuil[k] - 1, AC_sizes[k] / 4);
		Fourier::spectrum_by_autocorrelation_small_size(coords_per_mask[1][k], spectrum[1][k], masks_binar[k], seuil[k] - 1, AC_sizes[k] / 4);
		Fourier::spectrum_by_autocorrelation_small_size(coords_per_mask[2][k], spectrum[2][k], masks_binar[k], seuil[k] - 1, AC_sizes[k] / 4);
	}

	Synthesis(base_dir, filename_source,/*synthese_name,*//*N_input,*/S_input, pcas, spectrum, guidance_synthesis, coords_per_mask);
}

}
