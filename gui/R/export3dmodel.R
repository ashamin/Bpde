fx = c(x, rep(-1, length(H)-length(x)));
fy = c(y, rep(-1, length(H)-length(y)));
fH = c(H);
csv = data.frame(fx, fy, fH);
write.csv2(csv, file=csvFilePath);
