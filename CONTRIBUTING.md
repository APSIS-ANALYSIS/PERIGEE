# Contributing Guidelines

Thank you for your interest in contributing to the **PERIGEE** project! We welcome contributions from the community, whether through bug fixes, new features, or documentation improvements.

Please follow these guidelines to ensure that your contributions are as smooth and efficient as possible.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [How to Contribute](#how-to-contribute)
  - [Reporting Bugs](#reporting-bugs)
  - [Suggesting Enhancements](#suggesting-enhancements)
  - [Pull Requests](#pull-requests)
- [Coding Standards](#coding-standards)
- [Style Guidelines](#style-guidelines)
- [Testing](#testing)
- [Commit Messages](#commit-messages)
- [Licensing](#licensing)

## Code of Conduct

By participating in this project, you agree to be respectful and considerate to others.

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue in the [Issues](https://github.com/APSIS-ANALYSIS/PERIGEE/issues) section of this repository. Include as much detail as possible:
- Steps to reproduce the bug
- Expected behavior
- Actual behavior
- Any relevant logs or error messages
- Your environment (OS, compiler, versions, etc.)

### Suggesting Enhancements

We welcome suggestions for improvements! If you have an idea for a new feature or an improvement to an existing one, please open an issue and describe your idea. Be sure to include the context and the problem you're aiming to solve.

### Pull Requests

- **Fork the repository**: Create a personal fork of the project.
- **Clone your fork**: `git clone https://github.com/your-username/PERIGEE.git`
- **Create a new branch**: `git checkout -b feature/your-feature-name`
- **Make changes**: Work on your feature or bug fix.
- **Write tests**: If applicable, write tests for your changes.
- **Commit your changes**: Follow the [commit message guidelines](#commit-messages).
- **Push your changes**: `git push origin feature/your-feature-name`
- **Create a pull request**: Open a pull request from your branch to the `main` branch of this repository.

Please ensure your pull request follows these guidelines:
- **Be clear and concise**: Explain what changes you made and why.
- **Include tests**: Any new features or bug fixes should include tests.
- **Adhere to the code style**: Follow the [coding standards](#coding-standards) and [style guidelines](#style-guidelines).

## Coding Standards

- Use **C++11**.
- Follow consistent naming conventions and formatting.
- Use modern C++ features (e.g., `std::array`, `std::vector`, `auto`, smart pointers).
- Avoid `using namespace std;` in header files.
- Ensure that code is **readable** and **well-documented**.

## Style Guidelines

- **Indentation**: Use 2 spaces for indentation.
- **Braces**: Use braces `{}` for all control structures.
- **Line Length**: Limit lines to 80 characters for better readability.
- **Naming Conventions**: 
  - Classes should be named using PascalCase (e.g., `FEAElement`).
  - Functions and variables should use camelCase (e.g., `getElementType`).

## Testing

Before submitting your pull request, please ensure that:
- Your changes do not break existing functionality.
- All tests pass (if applicable).
- New features or bug fixes are accompanied by tests.

To run the tests locally:
1. Set up your environment according to the [installation instructions](INSTALL.md).
2. Run the test suite using the provided test framework.

## Commit Messages

Please follow these guidelines for your commit messages:
- **Keep the subject line short** (50 characters or less).
- **Use the imperative mood**: "Fix bug" instead of "Fixed bug."
- **Include a detailed description** for complex changes.
- **Reference issues**: If your commit addresses an issue, include the issue number, e.g., `Fixes #123`.
