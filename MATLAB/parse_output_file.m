function [errMat, dofMat, dtMat, setupMat, runMat, M_vec, p_vec] = parse_output_file(filename)

    fid = fopen(filename,'r');
    if fid < 0
        error('File could not be opened.');
    end

    % Read entire file into cell array of strings
    C = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    lines = C{1};

    % Helper function to extract a numeric vector from a line like:
    % p_vec: [1, 2, 3]
    function vec = extract_vec(str)
        vec_str = regexp(str, '\[(.*?)\]', 'tokens', 'once');
        nums = str2num(vec_str{1}); %#ok<ST2NM>
        vec = nums(:);
    end

    % Find indices for important sections
    idx_p = find(contains(lines, 'p_vec:'), 1);
    idx_M = find(contains(lines, 'M_vec:'), 1);
    idx_err = find(contains(lines, 'Error matrix'), 1);
    idx_dof = find(contains(lines, 'DOF matrix'), 1);
    idx_dt = find(contains(lines, 'dt matrix'), 1);
    idx_setup = find(contains(lines, 'Setup times'), 1);
    idx_run = find(contains(lines, 'Run times'), 1);

    % Extract size vectors
    p_vec = extract_vec(lines{idx_p});
    M_vec = extract_vec(lines{idx_M});
    nM = length(M_vec);
    nP = length(p_vec);

    % Extract a numeric matrix block (nM lines after header)
    function A = extract_matrix(start_idx)
        data = zeros(nM, nP);
        for i = 1:nM
            vals = sscanf(lines{start_idx + i}, '%f');
            data(i,:) = vals(1:nP).';
        end
        A = data;
    end

    % Parse matrices into double
    errMat   = extract_matrix(idx_err);
    dofMat   = extract_matrix(idx_dof);
    dtMat    = extract_matrix(idx_dt);
    setupMat = extract_matrix(idx_setup);
    runMat   = extract_matrix(idx_run);
end
