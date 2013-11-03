/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/tools/globals.h>

namespace votca { namespace tools {

	bool globals::verbose = false;
	std::string globals::url = "http://www.votca.org";
	std::string globals::email = "devs@votca.org";
        
        std::string globals::option_fmt(".TP\n\\fB%1%\\fR\n%2%\n");
        
        std::string globals::header_fmt(".TH \"%1%\" 1 \"\" \"Version: %2%\"\n\n");
        
        std::string globals::name_fmt(".SH NAME\n"
                           "\n.P\n"
                           "%1% \\- Part of the VOTCA package\n"
                           "\n.P\n"
                           "For more info please visit %2%\n\n"
        );
        
        std::string globals::authors_fmt ("\n.SH AUTHORS\n"
                            "\n.P\n"
                            "Written and maintained by the VOTCA Development Team <%1%>\n");
        
        std::string globals::copyright_fmt ("\n.SH COPYRIGHT\n"
                              "\n.P\n\n"
                              "Copyright 2009\\-2013 The VOTCA Development Team (%1%).\n"
                              "\n.P\n"
                              "Licensed under the Apache License, Version 2.0 (the \"License\") "
                              "you may not use this file except in compliance with the License. "
                              "You may obtain a copy of the License at"
                              "\n.P\n"
                              "http://www.apache.org/licenses/LICENSE\\-2.0\n"
                              "\n.P\n"
                              "Unless required by applicable law or agreed to in writing, software "
                              "distributed under the License is distributed on an \"AS IS\" BASIS, "
                              "WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. "
                              "See the License for the specific language governing permissions and "
                              "limitations under the License.");
               
        std::string globals::synopsis_fmt ("\n.SH SYNOPSIS\n"
                                "\n.P\n\\fB%1%\\fR [\\fIOPTION\\fR] [\\fIARGUMENT\\fR]\n");
        
        std::string globals::description_fmt ("\n.SH DESCRIPTION\n"
                                "\n.P\n%1%\n");
        
        std::string globals::options_fmt ("\n.SH OPTIONS\n");        

}}
