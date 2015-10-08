function marker_set_animation(S)

  fh = figure;
  ph = plot3(nan,nan,nan);

  for i=1:numel(S.time)
    [az,el] = view;
    if (mod(i+4,5) ~= 0)
      continue;
    end
    clf;
    tic
    plot3(S.Markers.CentreOfMass.X(i), ...
          S.Markers.CentreOfMass.Y(i), ...
          S.Markers.CentreOfMass.Z(i),'ro','MarkerSize',10), hold on

    % Right Leg
    plot3([S.Markers.RTOE.X(i), S.Markers.RHEE.X(i)], ...
          [S.Markers.RTOE.Y(i), S.Markers.RHEE.Y(i)], ...
          [S.Markers.RTOE.Z(i), S.Markers.RHEE.Z(i)], 'ko-', ...
          [S.Markers.RTOE.X(i), S.Markers.RANK.X(i)], ...
          [S.Markers.RTOE.Y(i), S.Markers.RANK.Y(i)], ...
          [S.Markers.RTOE.Z(i), S.Markers.RANK.Z(i)], 'ko-', ...
          [S.Markers.RHEE.X(i), S.Markers.RANK.X(i)], ...
          [S.Markers.RHEE.Y(i), S.Markers.RANK.Y(i)], ...
          [S.Markers.RHEE.Z(i), S.Markers.RANK.Z(i)], 'ko-', ...
          [S.Markers.RANK.X(i), S.Markers.RKNE.X(i)], ...
          [S.Markers.RANK.Y(i), S.Markers.RKNE.Y(i)], ...
          [S.Markers.RANK.Z(i), S.Markers.RKNE.Z(i)], 'ko-', ...
          [S.Markers.RTOE.X(i), S.Markers.RHEE.X(i)], ...
          [S.Markers.RTOE.Y(i), S.Markers.RHEE.Y(i)], ...
          [S.Markers.RTOE.Z(i), S.Markers.RHEE.Z(i)], 'ko-', ...
          [S.Markers.RPSI.X(i), S.Markers.RASI.X(i)], ...
          [S.Markers.RPSI.Y(i), S.Markers.RASI.Y(i)], ...
          [S.Markers.RPSI.Z(i), S.Markers.RASI.Z(i)], 'ko-', ...
          [S.Markers.RPSI.X(i), S.Markers.RKNE.X(i)], ...
          [S.Markers.RPSI.Y(i), S.Markers.RKNE.Y(i)], ...
          [S.Markers.RPSI.Z(i), S.Markers.RKNE.Z(i)], 'ko-', ...
          [S.Markers.RKNE.X(i), S.Markers.RASI.X(i)], ...
          [S.Markers.RKNE.Y(i), S.Markers.RASI.Y(i)], ...
          [S.Markers.RKNE.Z(i), S.Markers.RASI.Z(i)], 'ko-', ...
          [S.Markers.RKNE.X(i), S.Markers.RHEE.X(i)], ...
          [S.Markers.RKNE.Y(i), S.Markers.RHEE.Y(i)], ...
          [S.Markers.RKNE.Z(i), S.Markers.RHEE.Z(i)], 'ko-')

    % Left Leg
    plot3([S.Markers.LTOE.X(i), S.Markers.LHEE.X(i)], ...
          [S.Markers.LTOE.Y(i), S.Markers.LHEE.Y(i)], ...
          [S.Markers.LTOE.Z(i), S.Markers.LHEE.Z(i)], 'ko-', ...
          [S.Markers.LTOE.X(i), S.Markers.LANK.X(i)], ...
          [S.Markers.LTOE.Y(i), S.Markers.LANK.Y(i)], ...
          [S.Markers.LTOE.Z(i), S.Markers.LANK.Z(i)], 'ko-', ...
          [S.Markers.LHEE.X(i), S.Markers.LANK.X(i)], ...
          [S.Markers.LHEE.Y(i), S.Markers.LANK.Y(i)], ...
          [S.Markers.LHEE.Z(i), S.Markers.LANK.Z(i)], 'ko-', ...
          [S.Markers.LANK.X(i), S.Markers.LKNE.X(i)], ...
          [S.Markers.LANK.Y(i), S.Markers.LKNE.Y(i)], ...
          [S.Markers.LANK.Z(i), S.Markers.LKNE.Z(i)], 'ko-', ...
          [S.Markers.LTOE.X(i), S.Markers.LHEE.X(i)], ...
          [S.Markers.LTOE.Y(i), S.Markers.LHEE.Y(i)], ...
          [S.Markers.LTOE.Z(i), S.Markers.LHEE.Z(i)], 'ko-', ...
          [S.Markers.LPSI.X(i), S.Markers.LASI.X(i)], ...
          [S.Markers.LPSI.Y(i), S.Markers.LASI.Y(i)], ...
          [S.Markers.LPSI.Z(i), S.Markers.LASI.Z(i)], 'ko-', ...
          [S.Markers.LPSI.X(i), S.Markers.LKNE.X(i)], ...
          [S.Markers.LPSI.Y(i), S.Markers.LKNE.Y(i)], ...
          [S.Markers.LPSI.Z(i), S.Markers.LKNE.Z(i)], 'ko-', ...
          [S.Markers.LKNE.X(i), S.Markers.LASI.X(i)], ...
          [S.Markers.LKNE.Y(i), S.Markers.LASI.Y(i)], ...
          [S.Markers.LKNE.Z(i), S.Markers.LASI.Z(i)], 'ko-', ...
          [S.Markers.LKNE.X(i), S.Markers.LHEE.X(i)], ...
          [S.Markers.LKNE.Y(i), S.Markers.LHEE.Y(i)], ...
          [S.Markers.LKNE.Z(i), S.Markers.LHEE.Z(i)], 'ko-')

    % Right Arm
    plot3([S.Markers.RWRA.X(i), S.Markers.RWRB.X(i)], ...
          [S.Markers.RWRA.Y(i), S.Markers.RWRB.Y(i)], ...
          [S.Markers.RWRA.Z(i), S.Markers.RWRB.Z(i)], 'ko-', ...
          [S.Markers.RSHO.X(i), S.Markers.RELB.X(i)], ...
          [S.Markers.RSHO.Y(i), S.Markers.RELB.Y(i)], ...
          [S.Markers.RSHO.Z(i), S.Markers.RELB.Z(i)], 'ko-', ...
          [S.Markers.RWRA.X(i), S.Markers.RELB.X(i)], ...
          [S.Markers.RWRA.Y(i), S.Markers.RELB.Y(i)], ...
          [S.Markers.RWRA.Z(i), S.Markers.RELB.Z(i)], 'ko-', ...
          [S.Markers.RWRB.X(i), S.Markers.RELB.X(i)], ...
          [S.Markers.RWRB.Y(i), S.Markers.RELB.Y(i)], ...
          [S.Markers.RWRB.Z(i), S.Markers.RELB.Z(i)], 'ko-', ...
          [S.Markers.RFIN.X(i), S.Markers.RWRB.X(i)], ...
          [S.Markers.RFIN.Y(i), S.Markers.RWRB.Y(i)], ...
          [S.Markers.RFIN.Z(i), S.Markers.RWRB.Z(i)], 'ko-', ...
          [S.Markers.RFIN.X(i), S.Markers.RWRA.X(i)], ...
          [S.Markers.RFIN.Y(i), S.Markers.RWRA.Y(i)], ...
          [S.Markers.RFIN.Z(i), S.Markers.RWRA.Z(i)], 'ko-')


    % Left Arm
    plot3([S.Markers.LWRA.X(i), S.Markers.LWRB.X(i)], ...
          [S.Markers.LWRA.Y(i), S.Markers.LWRB.Y(i)], ...
          [S.Markers.LWRA.Z(i), S.Markers.LWRB.Z(i)], 'ko-', ...
          [S.Markers.LSHO.X(i), S.Markers.LELB.X(i)], ...
          [S.Markers.LSHO.Y(i), S.Markers.LELB.Y(i)], ...
          [S.Markers.LSHO.Z(i), S.Markers.LELB.Z(i)], 'ko-', ...
          [S.Markers.LWRA.X(i), S.Markers.LELB.X(i)], ...
          [S.Markers.LWRA.Y(i), S.Markers.LELB.Y(i)], ...
          [S.Markers.LWRA.Z(i), S.Markers.LELB.Z(i)], 'ko-', ...
          [S.Markers.LWRB.X(i), S.Markers.LELB.X(i)], ...
          [S.Markers.LWRB.Y(i), S.Markers.LELB.Y(i)], ...
          [S.Markers.LWRB.Z(i), S.Markers.LELB.Z(i)], 'ko-', ...
          [S.Markers.LFIN.X(i), S.Markers.LWRB.X(i)], ...
          [S.Markers.LFIN.Y(i), S.Markers.LWRB.Y(i)], ...
          [S.Markers.LFIN.Z(i), S.Markers.LWRB.Z(i)], 'ko-', ...
          [S.Markers.LFIN.X(i), S.Markers.LWRA.X(i)], ...
          [S.Markers.LFIN.Y(i), S.Markers.LWRA.Y(i)], ...
          [S.Markers.LFIN.Z(i), S.Markers.LWRA.Z(i)], 'ko-')

    % Head, back, chest
    plot3([S.Markers.RFHD.X(i), S.Markers.LFHD.X(i)], ...
          [S.Markers.RFHD.Y(i), S.Markers.LFHD.Y(i)], ...
          [S.Markers.RFHD.Z(i), S.Markers.LFHD.Z(i)], 'ko-', ...
          [S.Markers.RBHD.X(i), S.Markers.LBHD.X(i)], ...
          [S.Markers.RBHD.Y(i), S.Markers.LBHD.Y(i)], ...
          [S.Markers.RBHD.Z(i), S.Markers.LBHD.Z(i)], 'ko-', ...
          [S.Markers.RBHD.X(i), S.Markers.LFHD.X(i)], ...
          [S.Markers.RBHD.Y(i), S.Markers.LFHD.Y(i)], ...
          [S.Markers.RBHD.Z(i), S.Markers.LFHD.Z(i)], 'ko-', ...
          [S.Markers.RFHD.X(i), S.Markers.LBHD.X(i)], ...
          [S.Markers.RFHD.Y(i), S.Markers.LBHD.Y(i)], ...
          [S.Markers.RFHD.Z(i), S.Markers.LBHD.Z(i)], 'ko-', ...
          [S.Markers.RFHD.X(i), S.Markers.RBHD.X(i)], ...
          [S.Markers.RFHD.Y(i), S.Markers.RBHD.Y(i)], ...
          [S.Markers.RFHD.Z(i), S.Markers.RBHD.Z(i)], 'ko-', ...
          [S.Markers.LFHD.X(i), S.Markers.LBHD.X(i)], ...
          [S.Markers.LFHD.Y(i), S.Markers.LBHD.Y(i)], ...
          [S.Markers.LFHD.Z(i), S.Markers.LBHD.Z(i)], 'ko-', ...
          [S.Markers.C7.X(i), S.Markers.T10.X(i)], ...
          [S.Markers.C7.Y(i), S.Markers.T10.Y(i)], ...
          [S.Markers.C7.Z(i), S.Markers.T10.Z(i)], 'ko-', ...
          [S.Markers.CLAV.X(i), S.Markers.STRN.X(i)], ...
          [S.Markers.CLAV.Y(i), S.Markers.STRN.Y(i)], ...
          [S.Markers.CLAV.Z(i), S.Markers.STRN.Z(i)], 'ko-', ...
          [S.Markers.CLAV.X(i), S.Markers.C7.X(i)], ...
          [S.Markers.CLAV.Y(i), S.Markers.C7.Y(i)], ...
          [S.Markers.CLAV.Z(i), S.Markers.C7.Z(i)], 'ko-', ...
          [S.Markers.T10.X(i), S.Markers.STRN.X(i)], ...
          [S.Markers.T10.Y(i), S.Markers.STRN.Y(i)], ...
          [S.Markers.T10.Z(i), S.Markers.STRN.Z(i)], 'ko-')

    % Trunk Connect
    plot3([S.Markers.T10.X(i), S.Markers.LPSI.X(i)], ...
          [S.Markers.T10.Y(i), S.Markers.LPSI.Y(i)], ...
          [S.Markers.T10.Z(i), S.Markers.LPSI.Z(i)], 'ko-', ...
          [S.Markers.T10.X(i), S.Markers.LASI.X(i)], ...
          [S.Markers.T10.Y(i), S.Markers.LASI.Y(i)], ...
          [S.Markers.T10.Z(i), S.Markers.LASI.Z(i)], 'ko-', ...
          [S.Markers.T10.X(i), S.Markers.RPSI.X(i)], ...
          [S.Markers.T10.Y(i), S.Markers.RPSI.Y(i)], ...
          [S.Markers.T10.Z(i), S.Markers.RPSI.Z(i)], 'ko-', ...
          [S.Markers.T10.X(i), S.Markers.RASI.X(i)], ...
          [S.Markers.T10.Y(i), S.Markers.RASI.Y(i)], ...
          [S.Markers.T10.Z(i), S.Markers.RASI.Z(i)], 'ko-', ...
          [S.Markers.C7.X(i), S.Markers.RSHO.X(i)], ...
          [S.Markers.C7.Y(i), S.Markers.RSHO.Y(i)], ...
          [S.Markers.C7.Z(i), S.Markers.RSHO.Z(i)], 'ko-', ...
          [S.Markers.C7.X(i), S.Markers.LSHO.X(i)], ...
          [S.Markers.C7.Y(i), S.Markers.LSHO.Y(i)], ...
          [S.Markers.C7.Z(i), S.Markers.LSHO.Z(i)], 'ko-', ...
          [S.Markers.CLAV.X(i), S.Markers.RSHO.X(i)], ...
          [S.Markers.CLAV.Y(i), S.Markers.RSHO.Y(i)], ...
          [S.Markers.CLAV.Z(i), S.Markers.RSHO.Z(i)], 'ko-', ...
          [S.Markers.CLAV.X(i), S.Markers.LSHO.X(i)], ...
          [S.Markers.CLAV.Y(i), S.Markers.LSHO.Y(i)], ...
          [S.Markers.CLAV.Z(i), S.Markers.LSHO.Z(i)], 'ko-', ...
          [S.Markers.C7.X(i), S.Markers.RBHD.X(i)], ...
          [S.Markers.C7.Y(i), S.Markers.RBHD.Y(i)], ...
          [S.Markers.C7.Z(i), S.Markers.RBHD.Z(i)], 'ko-', ...
          [S.Markers.C7.X(i), S.Markers.LBHD.X(i)], ...
          [S.Markers.C7.Y(i), S.Markers.LBHD.Y(i)], ...
          [S.Markers.C7.Z(i), S.Markers.LBHD.Z(i)], 'ko-')

    axis equal
    xlim([-800,800])
    ylim([100,2000])
    zlim([0,1700])
    box on
    title([S.name, '   ', sprintf('%0.2f', S.time(i)), ' (s)'], ...
          'FontSize', 20,'interpreter','none')
%    view(az,el);
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    set(gca,'ztick',[])
    set(gca,'zticklabel',[])
    XL = get(gca, 'XLim');
    YL = get(gca, 'YLim');

    size([XL(1), XL(2), XL(2), XL(1)]);
    size([YL(1), YL(1), YL(2), YL(2)]);
    patch([XL(1), XL(2), XL(2), XL(1)], [YL(1), YL(1), YL(2), YL(2)], ...
          [0 0 0 0], 'FaceColor', [0.9 0.9 0.9]);
    drawnow
    if i==1
      pause
      [az,el] = view;
    end
    toc
  end
end

