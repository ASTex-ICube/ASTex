#include "texton_io.h"

namespace ASTex
{

bool import_texton_from_png(ASTex::ImageRGBd& texton, std::string file_path)
{
    return IO::loadu8_in_01(texton, file_path);
}

//Import a texton from a file generated by the original mathlab implem
bool import_texton(ASTex::ImageRGBd& texton , std::string file_path, ImageRGBd::PixelType &mean)
{

    std::ifstream infile;
    infile.open(file_path, std::ifstream::in);

    std::cout <<" OPENED FILE "<< std::endl;

    // First line: Should be TEXTON
    std::string lineString;
    std::getline(infile,lineString);
    if(lineString.compare("TEXTON")!=0)
    {
        std::cerr << "Not a Texton file" << std::endl;
        std::cout <<lineString<< std::endl;

        return false;
    }
    // Second line: Version
    std::getline(infile,lineString);
    int version = atoi(lineString.c_str());
    if (version !=1 )
    {
        std::cerr << "Version number not supported: only supported versions are 1" << std::endl;
        return false;
    }

    // 3rd line: order
    std::getline(infile,lineString);
    int order = atoi(lineString.c_str());
    // 4th line: mean
    for(int i=0; i<3; ++i)
    {
        std::getline(infile,lineString,' ');
        mean[i] = atof(lineString.c_str());
    }
    // 5th line: size
    std::getline(infile,lineString,' ');
    int M = atoi(lineString.c_str());


    std::getline(infile,lineString);
    int N = atoi(lineString.c_str());

    std::cout << "version: " << version << std::endl;
    std::cout << "order " << order << std::endl;
    std::cout << "mean: " << mean[0] << " " << mean[1] << " " << mean[2] << std::endl;
    std::cout << "size: " << M << " " << N << std::endl;

    texton.initItk(M,N,true);

    for (int y = 0; y < N; ++y)
    {
        for (int x = 0; x < M; ++x)
        {
            std::getline(infile,lineString,' ');
            for(int i=0; i<3; ++i)
                    texton.pixelAbsolute(x,y)[i] = atof(lineString.c_str());
        }
    }

    std::cout <<" Texton loaded: "<< file_path << std::endl;
    infile.close();

    return true;
}

bool import_texton(ASTex::ImageRGBd& texton , std::string file_path)
{
    ImageRGBd::PixelType mean;

    std::ifstream infile;
    infile.open(file_path, std::ifstream::in);

    std::cout <<"OPENED FILE"<< std::endl;

    // First line: Should be TEXTON
    std::string lineString;
    std::getline(infile,lineString);
    if(lineString.compare("TEXTON")!=0)
    {
        std::cerr << "Not a Texton file" << std::endl;
        std::cout <<lineString<< std::endl;

        return false;
    }
    // Second line: Version
    std::getline(infile,lineString);
    int version = atoi(lineString.c_str());
    if (version !=1 )
    {
        std::cerr << "Version number not supported: only supported versions are 1" << std::endl;
        return false;
    }

    // 3rd line: order
    std::getline(infile,lineString);
    int order = atoi(lineString.c_str());
    // 4th line: mean
    std::getline(infile,lineString,' ');
    std::cout << lineString.c_str() << std::endl;
    std::cout << std::to_string( std::stod(lineString.c_str())) << std::endl;
    mean[0] = std::stod(lineString.c_str());
    std::getline(infile,lineString,' ');
    mean[1] = std::stod(lineString.c_str());
    std::getline(infile,lineString);
    mean[2] = std::stod(lineString.c_str());
    // 5th line: size
    std::getline(infile,lineString, ' ');
    int M = atoi(lineString.c_str());

    std::getline(infile,lineString);
    int N = atoi(lineString.c_str());

    std::cout << "version: " << version << std::endl;
    std::cout << "order " << order << std::endl;
    std::cout << "mean: " << mean[0] << " " << mean[1] << " " << mean[2] << std::endl;
    std::cout << "size: " << M << " " << N << std::endl;

    texton.initItk(M,N,true);

    for (int y = 0; y < N; ++y)
    {
        for (int x = 0; x < M; ++x)
        {
            int i=0;
            for(; i<2; ++i)
            {
                std::getline(infile,lineString,' ');
                texton.pixelAbsolute(M-x-1, y)[i] = std::stod(lineString.c_str()) * std::sqrt(M*N) + mean[i];
            }
            std::getline(infile, lineString);
            texton.pixelAbsolute(M-x-1, y)[i] = std::stod(lineString.c_str()) * std::sqrt(M*N) + mean[i];
        }
    }

    std::cout <<" Texton loaded: "<< file_path << std::endl;
    infile.close();

    return true;
}

bool export_texton(ASTex::ImageRGBd& texton, std::string file_path,
                  int order,
                  ImageRGBd::PixelType &mean)
{
    std::ofstream outfile;
    outfile.open(file_path, std::ofstream::binary);

    //line 1 Version
    outfile << "TEXTON" << std::endl;

    //line 2 Version
    outfile << 1 << std::endl;

    //line 3 order
    outfile << order << std::endl;

    //line 4 Mean R G B
    outfile << mean[0] << ' ' << mean[1] << ' ' << mean[2] << std::endl;

    //line 5 M N
    outfile << texton.width() << ' ' << texton.height() << std::endl;

    for (int y = 0; y <  texton.height(); ++y)
    {
        for (int x = 0; x < texton.width(); ++x)
        {
            outfile << texton.pixelAbsolute(x,y)[0]
                    << ' '
                    << texton.pixelAbsolute(x,y)[1]
                    << ' '
                    << texton.pixelAbsolute(x,y)[2]
                    << std::endl;
        }
    }

    outfile.close();

    return true;
}

}

