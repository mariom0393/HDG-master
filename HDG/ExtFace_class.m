function infoFaces = ExtFace_class(infoFaces,X,T)
k=1;
q=1;
m=1;
for i = 1:length(infoFaces.extFaces)
    El_num = T(infoFaces.extFaces(i,1),4:6); %calls the in-edge nodes
    Nod_1 = X(El_num(1),:); 
    Nod_2 = X(El_num(2),:); 
    Nod_3 = X(El_num(3),:); 
    if Nod_1(2) == 0 || Nod_2(2) == 0 || Nod_3(2) == 0 % y=0
        infoFaces.extFaces_N(k,:) = infoFaces.extFaces(i,:);
        k = k+1;
    elseif Nod_1(1) == 1 || Nod_2(1) == 1 || Nod_3(1) == 1 % x=1
        infoFaces.extFaces_R(m,:) = infoFaces.extFaces(i,:);
        m=m+1;
    else
        infoFaces.extFaces_D(q,:) = infoFaces.extFaces(i,:);
        q = q+1;
    end
end
end

