# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 10:47:39 2021

@author: Dorian
"""
import numpy
from json import JSONEncoder,JSONDecoder
import json

class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)

class NumpyArrayDecoder(JSONDecoder):
    def default(self, obj):
        if isinstance(obj, list):
            return numpy.asarray(obj)
        return JSONEncoder.default(self, obj)


if __name__ == "__main__":
    
    numpyArrayOne = numpy.array([[11 ,22, 33], [44, 55, 66], [77, 88, 99]])
    numpyArrayTwo = numpy.array([[51, 61, 91], [121 ,118, 127]])
    
    # Serialization
    numpyData = {"arrayOne": numpyArrayOne, "arrayTwo": numpyArrayTwo}
    print("Original Data: \n")
    print(numpyData)
    print("\nSerialize NumPy array into JSON and write into a file")
    with open("numpyData.json", "w") as write_file:
        json.dump(numpyData, write_file, cls=NumpyArrayEncoder)
    print("Done writing serialized NumPy array into file")
    
    # Deserialization
    print("Started Reading JSON file")
    with open("numpyData.json", "r") as read_file:
        print("Converting JSON encoded data into Numpy array")
        decodedArray = json.load(read_file, cls=NumpyArrayDecoder)
        
    print("Re-Imported Data: \n")
    print(decodedArray)

    