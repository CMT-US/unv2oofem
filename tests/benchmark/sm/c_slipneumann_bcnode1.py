#   Copyright (c) 2017, Adam Sciegaj, Chalmers University of Technology
#
#   This library is free software; you can redistribute it and/or
#   modify it under the terms of the GNU Lesser General Public
#   License as published by the Free Software Foundation; either
#   version 2.1 of the License, or (at your option) any later version.
#
#   This library is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#
#   You should have received a copy of the GNU Lesser General Public
#   License along with this library; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

import numpy as np
import re
def giveUserDefBC(iCoords, iDofNum, iTime):
    #Read the strain history from the txt file produced by Matlab/Octave
    fil=open('C_slipNeumann_esgHistory.txt')
    inpFile=fil.readlines()
    fil.close

    esgRegex=re.compile(r'(\d+)\s+(-?\d\.\d{6}e.\d{2})\s+(-?\d\.\d{6}e.\d{2})\s+(-?\d\.\d{6}e.\d{2})\s+(-?\d\.\d{6}e.\d{2})\s+(-?\d\.\d{6}e.\d{2})\s+(-?\d\.\d{6}e.\d{2})\s+(-?\d\.\d{6}e.\d{2})\s+(-?\d\.\d{6}e.\d{2})\s+(-?\d\.\d{6}e.\d{2})\n?')

    strainList=[]
    slipList=[]
    slipGradList=[]
    for i in range(len(inpFile)):
        if esgRegex.findall(inpFile[i]) != []:
            time=int(esgRegex.findall(inpFile[i])[0][0])
            exx=float(esgRegex.findall(inpFile[i])[0][1])
            eyy=float(esgRegex.findall(inpFile[i])[0][2])
            gamxy=float(esgRegex.findall(inpFile[i])[0][3])
            sx=float(esgRegex.findall(inpFile[i])[0][4])
            sy=float(esgRegex.findall(inpFile[i])[0][5])
            gxx=float(esgRegex.findall(inpFile[i])[0][6])
            gyy=float(esgRegex.findall(inpFile[i])[0][7])
            gxy=float(esgRegex.findall(inpFile[i])[0][8])
            gyx=float(esgRegex.findall(inpFile[i])[0][9])
            strainList.append([time, exx, eyy, gamxy])
            slipList.append([time,sx,sy])
            slipGradList.append([time,gxx,gyy,gxy,gyx])

    strainHist=np.array(strainList)
    slipHist=np.array(slipList)
    slipGradHist=np.array(slipGradList)

    presVal=0.0
    iTime=int(iTime)
    strain = [ [strainHist[iTime-1][1], strainHist[iTime-1][3]/2], [strainHist[iTime-1][3]/2 , strainHist[iTime-1][2] ] ]
    slip = [ slipHist[iTime-1][1], slipHist[iTime-1][2]]
    slipGradient = [ [slipGradHist[iTime-1][1], slipGradHist[iTime-1][3]], [slipGradHist[iTime-1][4] , slipGradHist[iTime-1][2] ] ]

    if iDofNum == 1:
        presVal = 0
    elif iDofNum == 2:
        presVal = 0
    elif iDofNum == 31:
        presVal = -2*slipGradient[0][0] - slipGradient[0][1] - 3*slip[0]
    elif iDofNum == 32:
        presVal = -2*slipGradient[1][0] - slipGradient[1][1] - 3*slip[1]
    else:
        print('Warning: no rule for dof number: ', iDofNum)

    #	print 'presVal: ',presVal
    return presVal


