from lxml import etree
tree = etree.parse('dftgwbse.xml')

# etree.strip_attributes(tree, 'help')
# etree.strip_attributes(tree, 'default')
# etree.strip_attributes(tree, 'choices')
# etree.strip_attributes(tree, 'link')
# etree.strip_attributes(tree, 'unit')

def recursively_remove_empty(e):
    
    print('starting with ', e, str(e.text)) 
    
    c = list(recursively_remove_empty(c) for c in e.iterchildren())   
    print(c, all(c))
    try:
        if len(list(c)) == 0:
            if (e.text is None) or (e.text.strip() == ''):
                print('remove', e, e.text)
                e.getparent().remove(e)
                return True
            else:
                print('keep ',e, e.text)
                return False
        
        elif all(list(c)):
            print('remove all', e, list(c))
            e.getparent().remove(e)
            return True
    except:        
        print('issue with ',e)
        return False

recursively_remove_empty(tree.getroot())


# for elem in tree.xpath('//*[not(node())]'):
#      elem.getparent().remove(elem)



# for elem in tree.xpath('//*[* and not(*[*]) and not(*[not(*) and normalize-space()])]'):
#      elem.getparent().remove(elem)


# for elem in tree.xpath('//*[count(child::*) = 0]'):
    # elem.getparent().remove(elem)


# for elem in tree.xpath('//*[not(text())]'):
#     elem.getparent().remove(elem)




print(etree.dump(tree.getroot()))
tree.write('aaa.xml')