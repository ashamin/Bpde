svg(width=6,height=6,pointsize=10,filename=tfile);
library("fields");
image.plot(x, y, H);
contour(x, y, H, add = TRUE);
dev.off();
