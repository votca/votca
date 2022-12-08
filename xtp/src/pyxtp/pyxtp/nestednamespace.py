from types import new_class, SimpleNamespace
import xmltodict 
import os

class NestedNamespace(SimpleNamespace):
    """Extend SimpleNamespace to nested structure

    example:
    .. code-block:: python
       d = NestedNamespace({'a': 1, 'b':{'x':2, 'y':3}})
       d.a    # 1
       d.b.x  # 2

    In case of elements starting with a @ they are replaced
    by underscore, e.g. :

    .. code-block:: python
       d = NestedNamespace({'a': 1, 'b':{'@x':2, 'y':3}})
       d.a    # 1
       d.b._x  # 2
    """

    def __init__(self, dictionary, **kwargs):

        super().__init__(**kwargs)

        for key, value in dictionary.items():
            if isinstance(value, dict):
                self.__setattr__(key, NestedNamespace(value))
            else:
                if key.startswith('@'):
                    key = '_'+key[1:]
                self.__setattr__(key, value)
                
    @classmethod
    def __fromxml__(cls, xml_str: str, xml_attribs=True):
        """Creates a class instance from an xml file stored in a string

        Args:
            xml_str (str): content of an xml file represented as string

        Returns:
            NestedNamespace: namesapce mapping the xml file
            
        examle:
        .. code-block:: python
            from lxml import etree
            tree = etree.parse('test.xml')
            xml_str = etree.tostring(tree)
            data = NestedNamespace._fromxml(xml_str)
        """
        dict_data = xmltodict.parse(xml_str, xml_attribs=xml_attribs)
        return cls(dict_data)   
    
    def __todict__(self) -> dict:
        """Transform the nested namespace into a dictionary.
        
        Returns:
            dict: dictionary maps the namespace
        """
        
        def _recursive_namespace2dict(namespace_data: NestedNamespace, 
                        path: str='') -> dict:
            """Transform a nested namespace into a dictionary.

            Args:
                namespace_data (NestedNamespace): data in namespace form
                path (str, optional): basepath in the structure. Defaults to ''.
                remove_optional(bool, optional): remove the attr starting with '_' Default to True
            
            Returns:
                dict: dictionary maps the namespace
            """
            output = dict()
            for key, value in namespace_data.__dict__.items():
                if isinstance(value, NestedNamespace):
                    output.update(_recursive_namespace2dict(value, path=os.path.join(path,key)))
                else:
                    if not key.startswith('_'):
                        output[os.path.join(path,key)] = value
            return output

        return _recursive_namespace2dict(self)


def make_options(name: str, xml_dict: dict, set_default: bool=True):
    """Create a config object out of an XML dictionary

    Parameters
    ----------
    name : str
      Name for options type

    xml_dict : dict
      Dictionary generated from the XML configuration files

    set_default : bool (default: True)
      Whether to set default values in the config object

    Returns
    -------
    name_t

    """
    attrs, slots = {}, {}
    for key, value in xml_dict.items():
        if isinstance(value, dict):
            attrs[key] = make_options(key, value, set_default)
            doc = [v for k, v in value.items() if k.startswith("@")]
            slots[key] = "\n".join(doc)
        else:  # always a @key
            _key = key.replace("@", "_")
            attrs[_key] = value
            slots[_key] = ""
    cls = new_class(
        f"{name}_t", exec_body=lambda ns: ns.update({"__slots__": slots})
    )
    obj = cls()
    for k, v in attrs.items():
        setattr(obj, k, v)
        if not set_default:
            continue
        if k.startswith("_"):
            continue
        default = getattr(v, "_default", None)
        if default == "OPTIONAL":
            default = None
        children = [a for a in dir(v) if not a.startswith("_")]
        if len(children) == 0:
            setattr(obj, k, default)
    return obj
