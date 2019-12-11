function output_data(fdat,fmet,pstr,pnum,nstep,ctime,u,p,y)

% Write time data to binary file
fwrite(fdat,u,'double');
fwrite(fdat,p,'double');
fwrite(fdat,y,'double');

% Write time data to metadata file 
fprintf(fmet,'%6i %9.2e \n',nstep,ctime);

end