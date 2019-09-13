/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#include <expat.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stack>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <string>

#include <votca/tools/colors.h>
#include <votca/tools/property.h>
#include <votca/tools/propertyiomanipulator.h>
#include <votca/tools/tokenizer.h>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <unistd.h>

namespace votca {
namespace tools {
using namespace std;
// ostream modifier defines the output format, level, indentation
const int Property::IOindex = std::ios_base::xalloc();

const Property &Property::get(const string &key) const {
  Tokenizer tok(key, ".");
  Tokenizer::iterator n = tok.begin();
  if (n == tok.end()) return *this;

  const Property *p;
  map<string, int>::const_iterator iter;
  if (*n == "") {
    p = this;
  } else {
    iter = _map.find(*n);
    if (iter == _map.end())
      throw std::runtime_error("property not found: " + key);
    p = &_properties[((*iter).second)];
  }
  ++n;
  try {
    for (; n != tok.end(); ++n) {
      p = &p->get(*n);
    }
  } catch (std::runtime_error &err) {  // catch here to get full key in
                                       // exception
    throw std::runtime_error("property not found: " + key);
  }

  return *p;
}

Property &Property::get(const string &key) {
  return const_cast<Property &>(static_cast<const Property &>(*this).get(key));
}

Property &Property::getOradd(const std::string &key) {
  if (exists(key)) {
    return get(key);
  } else {
    return add(key, "");
  }
}

std::vector<const Property *> Property::Select(const string &filter) const {
  Tokenizer tok(filter, ".");
  std::vector<const Property *> selection;
  if (tok.begin() == tok.end()) return selection;
  selection.push_back(this);
  for (const auto &n : tok) {
    std::vector<const Property *> childs;
    for (const Property *p : selection) {
      for (const Property &s : p->_properties) {
        if (wildcmp(n.c_str(), s.name().c_str())) {
          childs.push_back(&s);
        }
      }
    }
    selection = childs;
  }
  return selection;
}

std::vector<Property *> Property::Select(const string &filter) {
  Tokenizer tok(filter, ".");
  std::vector<Property *> selection;
  if (tok.begin() == tok.end()) return selection;
  selection.push_back(this);
  for (const auto &n : tok) {
    std::vector<Property *> childs;
    for (Property *p : selection) {
      for (Property &s : p->_properties) {
        if (wildcmp(n.c_str(), s.name().c_str())) {
          childs.push_back(&s);
        }
      }
    }
    selection = childs;
  }
  return selection;
}

static void start_hndl(void *data, const char *el, const char **attr) {
  stack<Property *> *property_stack =
      (stack<Property *> *)XML_GetUserData((XML_Parser *)data);

  Property *cur = property_stack->top();
  Property &np = cur->add(el, "");

  for (int i = 0; attr[i]; i += 2) np.setAttribute(attr[i], attr[i + 1]);

  property_stack->push(&np);
}

static void end_hndl(void *data, const char *el) {
  stack<Property *> *property_stack =
      (stack<Property *> *)XML_GetUserData((XML_Parser *)data);
  property_stack->pop();
}

void char_hndl(void *data, const char *txt, int txtlen) {
  stack<Property *> *property_stack =
      (stack<Property *> *)XML_GetUserData((XML_Parser *)data);

  Property *cur = property_stack->top();
  cur->value().append(txt, txtlen);
}

void Property::LoadFromXML(string filename) {
  ifstream fl;
  fl.open(filename);
  if (!fl.is_open())
    throw std::ios_base::failure("Error on open xml file: " + filename);

  XML_Parser parser = XML_ParserCreate(NULL);
  if (!parser)
    throw std::runtime_error("Couldn't allocate memory for xml parser");

  XML_UseParserAsHandlerArg(parser);
  XML_SetElementHandler(parser, start_hndl, end_hndl);
  XML_SetCharacterDataHandler(parser, char_hndl);

  stack<Property *> pstack;
  pstack.push(this);

  XML_SetUserData(parser, (void *)&pstack);
  while (!fl.eof()) {
    string line;
    getline(fl, line);
    line = line + "\n";
    if (!XML_Parse(parser, line.c_str(), line.length(), fl.eof())) {
      throw std::ios_base::failure(
          filename + ": Parse error at line " +
          boost::lexical_cast<string>(XML_GetCurrentLineNumber(parser)) + "\n" +
          XML_ErrorString(XML_GetErrorCode(parser)));
    }
  }
  fl.close();
  XML_ParserFree(parser);
}

void PrintNodeTXT(std::ostream &out, const Property &p, const int start_level,
                  int level = 0, string prefix = "", string offset = "") {
  if ((p.value() != "") || p.HasChildren()) {
    if (level >= start_level) {
      if ((p.value()).find_first_not_of("\t\n ") != std::string::npos)
        out << offset << prefix << " = " << p.value() << endl;
    } else {
      prefix = "";
    }
  }

  for (const Property &prop : p) {
    level++;
    if (prefix == "") {
      PrintNodeTXT(out, prop, start_level, level, prefix + prop.name(), offset);
    } else {
      PrintNodeTXT(out, prop, start_level, level, prefix + "." + prop.name(),
                   offset);
    }
    level--;
  }
}

void PrintNodeXML(std::ostream &out, const Property &p,
                  PropertyIOManipulator *piom, int level = 0,
                  string offset = "") {
  Property::const_AttributeIterator ia;
  bool linebreak = true;
  bool has_value = false;

  const ColorSchemeBase *color = &DEFAULT_COLORS;

  string indent("");
  int start_level(0);

  if (piom) {
    start_level = piom->getLevel();
    indent = piom->getIndentation();
    color = piom->getColorScheme();
  }

  string cKey = color->Magenta();
  string cAttribute = color->Blue();
  string cAttributeValue = color->Green();
  string cReset = color->Reset();

  // print starting only from the start_level (the first node (level 0) can be
  // <> </>)
  if (level >= start_level) {
    // print the node name
    out << indent << offset << "<" << cKey << p.name() << cReset;
    // print the node attributes
    for (ia = p.firstAttribute(); ia != p.lastAttribute(); ++ia) {
      out << " " << cAttribute << ia->first << cReset << "=\""
          << cAttributeValue << ia->second << cReset << "\"";
    }
    // print node value if it is not empty
    has_value = ((p.value()).find_first_not_of("\t\n ") != std::string::npos);
    if (has_value || p.HasChildren()) {
      out << ">";
    } else {
      out << "/>" << std::endl;
    }
    if (has_value) {
      out << cAttributeValue << p.value() << cReset;
      linebreak = false;
    }

    // check if we need the end of the line or not
    if (!has_value && p.HasChildren()) out << std::endl;
    if (!has_value && !p.HasChildren()) linebreak = false;
  }

  // continue iteratively through the rest of the nodes
  for (const Property &prop : p) {
    level++;
    if (level > start_level) offset += "\t";
    PrintNodeXML(out, prop, piom, level, offset);
    if (level > start_level) offset.resize(offset.size() - 1);
    level--;
  }

  if (level >= start_level) {
    if (linebreak) {
      out << indent << offset << "</" << cKey << p.name() << cReset << ">"
          << std::endl;
    } else if (has_value) {
      out << "</" << cKey << p.name() << cReset << ">" << std::endl;
    }
  }
}

void PrintNodeTEX(std::ostream &out, const Property &p,
                  PropertyIOManipulator *piom, int level = 0,
                  string prefix = "") {

  int start_level = 0;
  if (piom) {
    start_level = piom->getLevel();
  }

  string head_name;
  string section("");  // reference of the description section in the manual
  string help("");
  // if this is the head node, print the header
  if (level == start_level) {

    string header_format(
        "\\subsection{%1%}\n"
        "\\label{%2%}\n%3%\n"
        "\\rowcolors{1}{invisiblegray}{white}\n"
        "{\\small\n "
        "\\begin{longtable}{m{3cm}|m{2cm}|m{1cm}|m{8cm}}\n"
        " option & default & unit & description\\\\\n\\hline\n");

    head_name = p.name();
    string label =
        "calc:" + head_name;  // reference of the xml file in the manual
    if (p.hasAttribute("section")) section = p.getAttribute<string>("section");
    if (p.hasAttribute("help")) help = p.getAttribute<string>("help");
    out << boost::format(header_format) % head_name % label % help;
    prefix = p.name();
  }

  if (level > start_level) {
    // if this node has children or a value or is not the first, start recursive
    // printing
    if ((p.value() != "" || p.HasChildren()) && level > -1) {
      string tex_name = boost::replace_all_copy(p.name(), "_", "\\_");
      string defaults("");  // default value if supplied
      if (p.hasAttribute("default"))
        defaults = p.getAttribute<string>("default");
      string unit("");  // unit, if supplied
      if (p.hasAttribute("unit")) unit = p.getAttribute<string>("unit");
      if (p.hasAttribute("help")) help = p.getAttribute<string>("help");

      string body_format(
          " \\hspace{%1%pt}\\hypertarget{%2%}{%3%} & %4% & %5% & %6% \\\\\n");

      out << boost::format(body_format) % int((level - start_level - 1) * 10) %
                 prefix % tex_name % defaults % unit % help;
    }
  }

  // continue iteratively through the rest of the nodes
  for (const Property &pp : p) {
    level++;
    if (prefix == "") {
      PrintNodeTEX(out, pp, piom, level, prefix);
    } else {
      PrintNodeTEX(out, pp, piom, level, prefix);
    }
    level--;
  }

  // if this is the head node, print the footer
  if (level == start_level) {
    string footer_format(
        "\\end{longtable}\n}\n"
        "\\noindent Return to the description of "
        "\\slink{%1%}{\\texttt{%2%}}.\n");

    out << boost::format(footer_format) % section % head_name;
  }
}

void PrintNodeHLP(std::ostream &out, const Property &p,
                  const int start_level = 0, int level = 0,
                  const string &prefix = "", const string &offset = "") {

  typedef Color<csRGB> ColorRGB;  // use the RGB palette
  ColorRGB RGB;                   // Instance of an RGB palette
  string fmt = "t|%1%%|15t|" + string(RGB.Blue()) + "%2%" +
               string(RGB.Green()) + "%|40t|%3%%|55t|" + string(RGB.Reset()) +
               "%4%\n";

  int leveloffset = level;
  string help("");
  // if this is the head node, print the header
  if (level == start_level) {
    string head_name = string(RGB.Magenta()) + p.name();
    if (p.hasAttribute("help")) {
      if (p.hasAttribute("help"))
        help = string(RGB.Red()) + p.getAttribute<string>("help");
      out << boost::format(" %1%: %|18t| %2%" + string(RGB.Reset()) + "\n") %
                 head_name % help;
    }
    leveloffset = 0;
    out << boost::format("%|3" + fmt) % "OPTION" % "DEFAULT" % "UNIT" %
               "DESCRIPTION";
  }

  if (level > start_level) {
    string ofmt = "%|" + boost::lexical_cast<string>(leveloffset) + fmt;
    string unit("");
    if (p.hasAttribute("unit")) unit = p.getAttribute<string>("unit");
    string defaults("");
    if (p.hasAttribute("default")) defaults = p.getAttribute<string>("default");
    if (p.hasAttribute("help")) help = p.getAttribute<string>("help");
    if (!unit.empty()) unit = "[" + unit + "]";
    if (!defaults.empty()) defaults = "(" + defaults + ")";

    string name = p.name();

    out << boost::format(ofmt) % name % defaults % unit % help;
  }

  for (const Property pp : p) {
    level++;
    if (prefix == "") {
      PrintNodeHLP(out, pp, start_level, level, pp.name(), offset);
    } else {
      PrintNodeHLP(out, pp, start_level, level, prefix + "." + pp.name(),
                   offset);
    }
    level--;
  }
}

std::ostream &operator<<(std::ostream &out, const Property &p) {
  if (!out.good()) return out;

  std::ostream::sentry sentry(out);

  if (sentry) {
    // get the property format object attached to the stream
    PropertyIOManipulator *pm =
        (PropertyIOManipulator *)out.pword(Property::getIOindex());

    string indentation("");
    int level = 0;

    PropertyIOManipulator::Type type = PropertyIOManipulator::XML;
    if (pm) {
      indentation = pm->getIndentation();
      level = pm->getLevel();
      type = pm->getType();
      // check if we > or >> to a file and remove color codes
      // if ( out.tellp() != -1 )  - not suitable for pipes
      if (!isatty(STDOUT_FILENO) || !isatty(STDERR_FILENO)) {
        pm->setColorScheme<csDefault>();
      }
    }

    switch (type) {
      default:
        PrintNodeTXT(out, p, level);
      case PropertyIOManipulator::XML:
        PrintNodeXML(out, p, pm);
        break;
      case PropertyIOManipulator::TXT:
        PrintNodeTXT(out, p, level, 0, "", indentation);
        break;
      case PropertyIOManipulator::TEX:
        PrintNodeTEX(out, p, pm);
        break;
      case PropertyIOManipulator::HLP:
        PrintNodeHLP(out, p, level, 0, "", indentation);
        break;
    }
  }

  return out;
};
}  // namespace tools
}  // namespace votca
