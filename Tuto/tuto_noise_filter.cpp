#include <ASTex/noise.h>
#include <ASTex/easy_io.h>

using namespace ASTex;

int main()
{
    using T = double;
    ImageGray<T> im(512,512);
    Noise<T> noise;

    im.parallel_for_all_pixels([&](T &p,int i,int j)
    {
        T born_sup(2);
        T center = born_sup / T(2);
        T x = T(i) / T(im.width()) * born_sup - center;     // x between [-center,center]
        T y = T(j) / T(im.height()) * born_sup - center;    // y between [-center,center]
        T r = noise.basic2D(x,y);                           // noise value at (x,y)

        p = r;
    });

    IO::save01_in_u8(im,TEMPO_PATH + "noise.png");

    return 0;
}
