% NNSearch2DFEX query a GL-tree for the nearest neighbor(NN)
%
% SYNTAX
% [ NNG ] = NNSearch2DFEX( rp, qp );       
% [ NNG, Dist ] = NNSearch2DFEX( rp, qp );  
%
% INPUT PARAMETERS
%       rp: [Nrx2] double vectors coordinates of reference points
%       qp: [Nqx2] double vectors coordinates of query points
%
% OUTPUT PARAMETERS
% 	NNG:  [Nqx1] column vector, each rows contains the NN index
%         in the reference points. So row one is the NN to first query
%         point.
%   Dist: [Nqx1] array, Facultative output, each rows contains the
%                   distance values of the found NN.
% 
% GENERAL INFORMATIONS
%         -This function is faster if all query points are given once
%         instead of looping and pass one point each loop.
%
% Remember the rule: for every query points, find its closest reference
%                    point
%
% For question, suggestion, bug reports
% giaccariluigi@msn.com
% 
% Visit: <a href="http://www.advancedmcode.org/gltree.html"> The GLTree Web Page</a>
% 
% Author :     Luigi Giaccari & Zexi Liu
% Created :    10/10/2008
% Last Update: 07/01/2012


