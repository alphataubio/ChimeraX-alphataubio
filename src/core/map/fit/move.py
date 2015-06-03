# -----------------------------------------------------------------------------
# Remember model and atom positions relative to a given model.
#
class Position_History:

    def __init__(self):

        self.cur_position = -1
        self.positions = []
        self.max_positions = 100

    def record_position(self, models, atoms, base_model):

        plist = self.positions
        c = self.cur_position
        if c < len(plist)-1:
            del plist[c+1:]      # Erase redo states.
        ps = position_state(models, atoms, base_model)
        self.positions.append(ps)
        if len(plist) > self.max_positions:
            del plist[0]
        else:
            self.cur_position += 1

    def undo(self):

        while self.can_undo():
            ps = self.positions[self.cur_position]
            self.cur_position -= 1
            if restore_position(ps):
                return True
        return False

    def can_undo(self):
        return self.cur_position >= 0

    def redo(self):

        while self.can_redo():
            self.cur_position += 1
            ps = self.positions[self.cur_position]
            if restore_position(ps):
                return True
        return False

    def can_redo(self):
        return self.cur_position+1 < len(self.positions)

# -----------------------------------------------------------------------------
#
def position_state(models, atoms, base_model):

  btfinv = base_model.position.inverse()
  mset = set(models)
  mset.update(tuple(atoms.unique_molecules))
  model_transforms = []
  for m in mset:
      model_transforms.append((m, btfinv * m.position))
  atom_positions = (atoms, atoms.coords.copy())
  return (base_model, model_transforms, atom_positions)

# -----------------------------------------------------------------------------
#
def restore_position(pstate, angle_tolerance = 1e-5, shift_tolerance = 1e-5):

    base_model, model_transforms, atom_positions = pstate
    if base_model.was_deleted:
        return False
    changed = False
    for m, mtf in model_transforms:
        if m.was_deleted:
            continue
        tf = base_model.position * mtf
        if not tf.same(m.position, angle_tolerance, shift_tolerance):
            changed = True
        m.set_place(tf)

    atoms, xyz = apos
    if not (atoms.coords == xyz).all():
        changed = True
        atoms.coords = xyz
    return changed

# -----------------------------------------------------------------------------
#
def move_models_and_atoms(tf, models, atoms, move_whole_molecules, base_model):

    if move_whole_molecules:
        models = list(models) + list(atoms.unique_molecules)
        from ...molecule import Atoms
        atoms = Atoms()
    global position_history
    position_history.record_position(models, atoms, base_model)
    for m in models:
        m.set_position(tf * m.position)
    if len(atoms) > 0:
        atoms.coords = tf * atoms.coords
    position_history.record_position(models, atoms, base_model)

# -----------------------------------------------------------------------------
#
position_history = Position_History()
