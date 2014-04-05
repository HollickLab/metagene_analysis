#!/usr/bin/python
"""metageneMethods represent commonly called functions in metagene analysis

Requires:
    python 2 (https://www.python.org/downloads/)

Joy-El R.B. Talbot Copyright (c) 2014

The MIT License (MIT)

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files 
(the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from MetageneError import MetageneError

def confirm_integer(value, descriptor, minimum=None, maximum=None):
    """Confirm integer statust and potentially minimum or maximum bound.  Return boolean.
    
    Keyword Arguments:
    value -- value (any type) to test
    descriptor -- descriptor of value to add to error message
    minimum -- minimum value allowed (default None)
    maximum -- maximum value allowed (default None)
    """
    try:
        if value != int(value):
            raise MetageneError("{} is not an integer".format(descriptor))
    except ValueError:
        raise MetageneError("{} is not an integer".format(descriptor))

    above_minimum = True
    below_maximum = True
    if minimum is not None and value < minimum:
        above_minimum = False
    if maximum is not None and value > maximum:
        below_maximum = False

    if above_minimum and below_maximum:
        return True
    else:
        if not above_minimum and not below_maximum:
            raise MetageneError("{} is outside of boundaries: {}-{}".format(
                descriptor, minimum, maximum))
        elif not above_minimum:
            raise MetageneError("{} is less than minimum: {}".format(
                descriptor, minimum))
        else:
            raise MetageneError("{} is greater than maximum: {}".format(
                descriptor, maximum))
    # end of confirm_integer function
    
    
    
    
    
