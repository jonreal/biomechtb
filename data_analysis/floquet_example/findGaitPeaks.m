function [locs] = findGaitPeaks(x,threshold,temporal_thres)

  % Find Peaks and thresold
  [pks,locs] = findpeaks(x);
  indx = find(pks>threshold(1) & pks<threshold(2));
  locs = locs(indx);

  numel(locs)
  % Remove peaks close together
  indx = find(locs(2:end) - locs(1:end-1) < temporal_thres);
  locs(indx) = [];

end
