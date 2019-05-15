function hMax = computeMeshSizeTri2D(T,X,nOfElements)

hMax = 0;

for iElem=1:nOfElements
    iVert = T(iElem,:);
    xVert = X(iVert,2);
    v1 = [1 1 3];
    v2 = [2 3 2];
    
    hElem = sqrt(sum((xVert(v1,:) - xVert(v2,:)).^2,2));
    
    hMax = max(hMax,max(hElem));
    
end