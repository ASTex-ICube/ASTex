#include <string>

namespace ASTex
{
/**
 * @brief KMnoise
 * @param filename_source
 * @param base_dir
 * @param nb_clusters
 * @param fft_size
 */
void KMnoise(const std::string& filename_source, const std::string& base_dir, const int nb_clusters, const int fft_size);

}
