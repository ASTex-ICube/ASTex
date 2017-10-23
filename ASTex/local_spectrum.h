#ifndef __ASTEX_LOCAL_SPECTRUM__
#define __ASTEX_LOCAL_SPECTRUM__

#include <vector>
//#include <ASTex/image_rgb.h>
#include <ASTex/image_gray.h>
#include <ASTex/image_rgb.h>
#include <ASTex/image_spectral.h>


namespace ASTex
{

namespace Fourier
{

/**
* @brief computation and storage of local spectra
*/
class ASTEX_API LocalSpectrum
{
	public :
		inline LocalSpectrum () {}

		LocalSpectrum (const ImageGrayd& input, int windowsize = 32, char windowtype = 'h', std::string name = "lsp");

		LocalSpectrum (const ImageGrayd& input, int windowsize = 32, int welch_windowsize=64, int welch_step = 8);

		void welch_post_loading(int welch_windowsize=64, int welch_step = 8);
		 
		void distance_map_to_pixel_linear_weights (ImageGrayd& output, int i, int j);

		void distance_map_to_pixel_uniform_weights (ImageGrayd& output, int i, int j);

		void distance_map_to_spectum_linear_weights(ImageGrayd& output,const ImageSpectrald& sp_ref);

		void distance_map_to_spectum_uniform_weights(ImageGrayd& output,const ImageSpectrald& sp_ref);

		void compute_local_guidance_map (int filtre_size = 32);
		
		void mosaic (ImageSpectrald& output);

		void distance_mosaic (ImageGrayd& output, int step);
		
		char get_windowtype ();

		char get_windowsize ();

		std::string get_name ();

		static double g(double x, double y, double sigma);

		static double g_non_norm(double x, double y, double sigma);

		static double hamming(double x, double y, double window_size);

		inline ImageGrayd& local_guidance(int i, int j)
		{ return m_local_guidance[i + j*m_input_width];}

		inline const ImageGrayd& local_guidance(int i, int j) const
		{ return m_local_guidance[i + j*m_input_width];}


		inline ImageSpectrald& spectrum(int i, int j)
		{ return m_spectrum[i + j*m_input_width]; }

		inline const ImageSpectrald& spectrum(int i, int j) const
		{ return m_spectrum[i + j*m_input_width]; }

		inline int get_used_filter_size()
		{
			if (m_local_guidance.empty())
				return -1;
			return m_local_guidance[0].width();
		}

	private:
		int m_input_width;
		int m_input_height;
		int m_windowsize;
		int m_left_bnd;
		int m_right_bnd;
		char m_windowtype;
		std::string m_name;

		std::vector<ImageSpectrald> m_spectrum;
		std::vector<ImageGrayd> m_local_guidance;

};




}
}


#endif 
