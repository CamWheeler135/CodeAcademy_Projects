// Header file for declaring classes and functions for Tic Tac Toe game.
#include <string>

class Board {
public:
    // Member variables.
    char board[9];
    bool complete;
    int rounds;

    // Constructor declaration.
    Board();

    // Collect the user input.
    int collect_input(char player);

    // Check the move is legal.
    bool check_legal_move(char *board, int move);

    // Apply the input to the board / update the board with the user input.
    void update_board(char *board, int move, char player);

    // Display the board.
    void display_board(char *board);

    // Check for a winner.
    bool check_for_winner(char *board);
};