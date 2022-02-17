clear all;
close all;
clc;

bag = rosbag('objBagFiles/fourthTest.bag');
selection = select(bag,'Topic','/tf');
msg = readMessages(selection,'DataFormat','struct');

idx = 1;
for ii = 1:size(msg,1)
    if size(struct2table(msg{ii}.Transforms), 1) == 0
        continue;
    end

    objPose(idx, 1) = msg{ii}.Transforms.Transform.Translation.X;
    objPose(idx, 2) = msg{ii}.Transforms.Transform.Translation.Y;
    objPose(idx, 3) = msg{ii}.Transforms.Transform.Translation.Z;
    idx = idx + 1;
end

figure();
plot3(objPose(:,1), objPose(:,2), objPose(:,3), 'LineWidth', 1.5);
