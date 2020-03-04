#include "histogram.h"
#include "gaussian_transfer.h"
#include <unistd.h>
#include <cstdio>
#include <getopt.h>
#include <ASTex/easy_io.h>

using namespace ASTex;

static void usage(const char *msg = nullptr, ...) {
    if (msg) {
        va_list args;
        va_start(args, msg);
        fprintf(stderr, "generator: ");
        vfprintf(stderr, msg, args);
        fprintf(stderr, "\n");
    }
    fprintf(stderr, R"(usage: generator <command> [options] <filenames...>
    commands: histo, gaussianize

    histo options:
        -o              name of the output file (default histo).
        --outfile       name of the output file (default histo).
        -n              nb_bins (default 256).
        --nb_bins       nb_bins (default 256).
        -y              range max of y of the hidtogram (default 1.0).
        --y_range_max   range max of y of the hidtogram (default 1.0).
        -g              to get the gray histogram (default RGB).
        --gray          to get the gray histogram (default RGB).
        --16            to load an 16 bits image (default load 8 bits image).

    gaussianize options:
        -o              name of the output file (default histo).
        --outfile       name of the output file (default histo).
        -g              to transform a gray image (default transform RGB image).
        --gray          to transform a gray image (default transform RGB image).
        --16            to load and save an 16 bits image (default load and save 8 bits image).
)");
    exit(EXIT_FAILURE);
}

int histo(int argc, char **argv)
{
    if(argc < 2) usage(argv[0]);

    bool gray = false;
    bool _16bits = false;
    int nb_bins = 256;
    std::string filename("histo");
    float ymax = 1.0;

    int option_index = 0;
    static struct option long_options[] = {
        {"outfile",     required_argument,  nullptr,    'o'},
        {"nb_bins",     required_argument,  nullptr,    'n'},
        {"y_range_max", required_argument,  nullptr,    'y'},
        {"gray" ,       no_argument,        nullptr,    'g' },
        {"16",          no_argument,  nullptr,    16}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "o:n:y:g",long_options,&option_index)) != -1) {
        switch (opt) {
            case 16:
                _16bits = true;
                break;
            case 'g':
                gray = true;
                break;
            case 'n':
                nb_bins = std::atoi(optarg);
                break;
            case 'o':
                filename = optarg;
                break;
            case 'y':
                ymax = float(std::atof(optarg));
                break;
            default:
                usage(argv[0]);

        }

    }

    if(gray){
        Histogram<ImageGrayd> h;
        ImageGrayd img;
        if(_16bits)
            IO::loadu16_in_01(img, argv[optind]);
        else
            IO::loadu8_in_01(img, argv[optind]);
        h.computeHisto(img, nb_bins);
        h.exportHisto(filename, ymax);
    }
    else {
        Histogram<ImageRGBd> h;
        ImageRGBd img;
        if(_16bits)
            IO::loadu16_in_01(img, argv[optind]);
        else
            IO::loadu8_in_01(img, argv[optind]);
        h.computeHisto(img, nb_bins);
        h.exportHisto(filename, ymax);
    }

    return EXIT_SUCCESS;
}

int gaussianize(int argc, char **argv)
{
    if(argc < 2) usage(argv[0]);

    bool gray = false;
    bool _16bits = false;
    std::string filename("out.png");

    int option_index = 0;
    static struct option long_options[] = {
        {"outfile",     required_argument,  nullptr,    'o'},
        {"gray" ,       no_argument,        nullptr,    'g' },
        {"16",          no_argument,  nullptr,    16}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "o:g",long_options,&option_index)) != -1) {
        switch (opt) {
            case 16:
                _16bits = true;
                break;
            case 'g':
                gray = true;
                break;
            case 'o':
                filename = optarg;
                break;
            default:
                usage(argv[0]);

        }

    }

    if(gray){
        ImageGrayf img;
        if(_16bits)
            IO::loadu16_in_01(img, argv[optind]);
        else
            IO::loadu8_in_01(img, argv[optind]);

        ImageGrayf imgT(img.width(), img.height());
        Gaussian_transfer::ComputeTinput(img,imgT);

        if(_16bits)
            IO::save01_in_u16(imgT, filename);
        else
            IO::save01_in_u8(imgT, filename);
    }
    else {
        ImageRGBf img;
        if(_16bits)
            IO::loadu16_in_01(img, argv[optind]);
        else
            IO::loadu8_in_01(img, argv[optind]);

        ImageRGBf imgT(img.width(), img.height());
        Gaussian_transfer::ComputeTinput(img,imgT);

        if(_16bits)
            IO::save01_in_u16(imgT, filename);
        else
            IO::save01_in_u8(imgT, filename);
    }

    return EXIT_SUCCESS;

}

int main(int argc, char **argv){
    if(argc < 2) usage();

    if(!strcmp(argv[1], "histo"))
        return histo(argc-1,argv+1);
    else if (!strcmp(argv[1], "gaussianize"))
        return gaussianize(argc-1, argv+1);
    else usage("unknown command \"%s\"", argv[1]);

    return EXIT_SUCCESS;
}
