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



#include <string>
#include <iostream>
#include <fstream>
#include <set>


inline bool file_exist(const std::string& name) 
{
	std::ifstream f(name.c_str());
	return f.good();
}

inline bool is_dll(const std::string& mot)
{
	return (mot.rfind(".dll") != std::string::npos);
}

inline void backslashize(std::string& mot)
{
	for (std::size_t i=0;i< mot.length(); ++i)
		if (mot[i] == '/')
			mot[i] = '\\';
}

inline std::string extract_path(std::string& filename)
{
	std::size_t p = filename.rfind('\\');
	if (p == std::string::npos)
		p = filename.rfind('/');

	if (p == std::string::npos)
		return ".";

	return filename.substr(0, p);
}

int main(int argc, char **argv)
{
	char buffer[512];

	if (argc != 2)
	{
		std::cout << "Drop exec over me to copy necessary dlls in its folder" << std::endl;
		std::cin.getline(buffer, 512);
		return EXIT_FAILURE;
	}


	std::string inst = extract_path(std::string(argv[0]));
	std::string dumpbin_exe= "\"" +DUMPBIN + "\"";
	std::string dumpbin = dumpbin_exe + " /DEPENDENTS ";
	std::string tempo_dll_list_file = TEMP_DIR+"/tempo_astex_dll_deploy.tmp";
	std::string exec(argv[1]);
	std::string path_exec = extract_path(exec);

	std::string cmd = dumpbin + exec + " > " + tempo_dll_list_file;
	system(cmd.c_str());

	// set of name of used dll
	std::set<std::string> dlls;

	// first pass dlls used by exec
	std::ifstream f(tempo_dll_list_file.c_str());
	for (int i = 0; i<9; ++i)
		f.getline(buffer, 512);
	std::string mot;
	while (!f.eof())
	{
		f >> mot;
		if (is_dll(mot) && file_exist(inst+"\\"+mot))
			dlls.insert(mot);
	}
	f.close();

	// other passes dlls used by dlls
	bool nfinished = true;
	while (nfinished)
	{
		nfinished = false;
		for (const auto& dll : dlls)
		{
			std::string cmd = dumpbin + inst+"/"+dll + "  > " + tempo_dll_list_file;
			system(cmd.c_str());
			std::ifstream f(tempo_dll_list_file.c_str());
			for (int i=0;i<9;++i) 
				f.getline(buffer,512);
			std::string mot;
			while (!f.eof())
			{
				f >> mot;
				if (is_dll(mot))
					if (dlls.find(mot) == dlls.end())
						if (file_exist(inst + "/" + mot))
						{
							dlls.insert(mot);
							nfinished = true;
						}
			}
			f.close();
		}
	}

	remove(tempo_dll_list_file.c_str());

	// copy the dll
	backslashize(path_exec);
	backslashize(inst);
	for (const auto& dll : dlls)
	{
		std::cout << "copying " << dll << std::endl;
		std::string copy = std::string("copy /Y ") + inst + "\\" + dll + " " + path_exec +" >NUL";
		system(copy.c_str());
	}
	std::cout << dlls.size() << " dlls copied" << std::endl;

	std::cout << "hit a key to finish " << std::endl;
	std::cin.getline(buffer,512);

    return EXIT_SUCCESS;
}
