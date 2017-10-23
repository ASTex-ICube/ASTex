
#include <fstream>
#include <cmath>
#include <algorithm>

#include <ASTex/local_spectrum.h>
#include <ASTex/dll.h>


namespace ASTex
{


void ASTEX_API bilateral_filter(const ImageGrayd& input, ImageGrayd& output, int filtre_size, double id, double cd);

void ASTEX_API bilateral_filter(const ImageRGBd& input, ImageRGBd& output, int filtre_size, double id, double cd);



namespace Fourier
{

/**
 * @brief Filtre bilateral avec descripteur de frequence
 */
void ASTEX_API frequency_bilateral_filter(const ImageGrayd& input, ImageGrayd& output, int filtre_size, LocalSpectrum& lsp , double id, double cd);

void ASTEX_API frequency_bilateral_filter(const ImageRGBd& input_color, /*const ImageGrayd& input,*/ ImageRGBd& output, int filtre_size, LocalSpectrum& lsp , double id, double cd);

void ASTEX_API frequency_joint_bilateral_filter(const ImageGrayd& input, ImageGrayd& output,const std::vector<ImageGrayd>& guidance_maps, int filtre_size, double id, double cd);


} //namespace Fourier


} //namespace ASTex
