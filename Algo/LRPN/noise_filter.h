#include <string>

namespace ASTex
{

/**
 * @brief RPnoise
 * @param inputfile
 * @param path_out (must finish with / and exist)
 * @param size
 * @param size_fft
 * @param sig_freq
 */
void RPnoise(const std::string& inputfile, const std::string& path_out,int size, int size_fft, float sig_freq, int nb_step_filt);


}
