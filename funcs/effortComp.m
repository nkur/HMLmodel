function [drivEff, exploEff] = effortComp(CHat, deltaU)

    projMatrix = CHat' * inv(CHat * CHat') * CHat;
    drivEff = norm(projMatrix * deltaU, 2);     % Efforts expended towards row space of C_hat

%     exploEff = norm(deltaU - projCu, 2);

    nullChat = null(CHat);     projNullChat = nullChat * inv(nullChat' * nullChat) * nullChat';
    exploEff = norm(projNullChat * deltaU, 2);     % Efforts expended towards null space of C_hat

end