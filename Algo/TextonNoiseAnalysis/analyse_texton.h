
#include <ASTex/easy_io.h>


namespace ASTex
{

int import_texton(const std::string& file_path, ImageRGBd& texton );

/**
 * @brief decomposition periodic + smooth
 * @param input input grayd image
 * @return a pair of image (periodic,smooth)
 */
std::pair<ImageGrayd,ImageGrayd> decompo(const ImageGrayd& input);

/**
 * @brief corrige les images ?
 * @param in input
 * @param per periodic (clamped to [0,1])
 * @param smo smooth (shifted +0.5 & clamped to [0,1])
 */
void corrige(const ImageGrayd& in, ImageGrayd& per, ImageGrayd& smo);

/**
 * @brief analyse_texton
 * @param input_color
 * @param images_alpha_visu
 * @param images_alpha_visu_bounded
 */
void analyse_texton( const ASTex::ImageRGBd& input_color,
					 std::vector<ASTex::ImageRGBd>& images_alpha_visu,
					 std::vector<ASTex::ImageRGBd>& images_alpha_visu_bounded);


}
