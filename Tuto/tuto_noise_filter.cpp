#include <ASTex/image_gray.h>
#include <ASTex/noise.h>
#include <ASTex/easy_io.h>

using namespace ASTex;

int main()
{
    using T = double;
    ImageGray<T> im(512,512);
    Noise<T> noise;

    im.parallel_for_all_pixels([&](T &p,int x,int y)
    {
        T r = noise.basic2D(T(x),T(y)); // noise value at (x,y)

        p = r;
    });

    IO::save01_in_u8(im,TEMPO_PATH + "noise.png");

    return 0;
}
