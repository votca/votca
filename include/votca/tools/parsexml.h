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

#ifndef VOTCA_TOOLS_PARSEXML_H
#define VOTCA_TOOLS_PARSEXML_H

#include <list>
#include <map>
#include <stack>
#include <string>

namespace votca {
namespace tools {

/**
    \brief XML SAX parser (wrapper for expat)

    This class is a wrapper for the expat SAX interface. Parsing the xml file
    is done via callback functions (Element handlers).

    The class implements contains a stupid implementation for functors to
    handle callbacks to member functions but there is lots of room for
   improvement. So far only callbacks for start element handlers is implemented.
   The extension to EndElement handler (to signal if an element is close) is
   similar and straight forward. So far it was not needed and was therefore not
   done.

  */
class ParseXML {
 public:
  /// constructor
  ParseXML() {}
  /// destructor
  ~ParseXML() {}

  /**
   * \brief open an XML file and start parsing it
   * @param _filename file to read
   *
   * This functions opens a file and parses the nodes of an XML file.
   * Make sure to set the corresponding element handler before to redirect
   * element handling.
   */
  void Open(const std::string &filename);

  /**
   * \brief Set handler for next element (only member functions possible)
   * @param object instance of class which for callback
   * @param fkt member function for callback
   *
   * This function always has to be called after processing an element to
   * say what is coming next. Optionally call IgnoreElement
   */
  template <typename T>
  void NextHandler(T *object,
                   void (T::*fkt)(const std::string &,
                                  std::map<std::string, std::string> &));

  /**
   * \brief Ignore the content of this elements and all of its childs
   */
  void IgnoreChilds();

 private:
  // virtual void ParseRoot(const string &el, map<string, string> &attr);
  void ParseIgnore(const std::string &el,
                   std::map<std::string, std::string> &attr);

  /// end element callback for xml parser
  void StartElemHndl(const std::string &el,
                     std::map<std::string, std::string> &attr);
  /// end element callback for xml parser
  void EndElemHndl(const std::string &el);

  class Functor {
   public:
    Functor() {}
    virtual void operator()(const std::string &,
                            std::map<std::string, std::string> &) = 0;
    virtual ~Functor(){};
  };

  template <typename T>
  class FunctorMember : public Functor {
   public:
    typedef void (T::*fkt_t)(const std::string &,
                             std::map<std::string, std::string> &);

    FunctorMember(T *cls, fkt_t fkt) : _cls(cls), _fkt(fkt) {}

    void operator()(const std::string &el,
                    std::map<std::string, std::string> &attr) {
      (_cls->*_fkt)(el, attr);
    }

   private:
    T *_cls;
    fkt_t _fkt;
  };

  std::stack<Functor *> _stack_handler;
  Functor *_handler;

  friend void start_hndl(void *data, const char *el, const char **attr);
  friend void end_hndl(void *data, const char *el);
};

inline void ParseXML::IgnoreChilds() {
  NextHandler(this, &ParseXML::ParseIgnore);
}

template <typename T>
inline void ParseXML::NextHandler(
    T *object,
    void (T::*fkt)(const std::string &, std::map<std::string, std::string> &)) {
  _handler = dynamic_cast<Functor *>(new FunctorMember<T>(object, fkt));
  _stack_handler.push(_handler);
}

}  // namespace tools
}  // namespace votca

#endif
